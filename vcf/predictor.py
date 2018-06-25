from collections import defaultdict, Counter
from typing import Any, Dict, List, Set
from itertools import combinations

from vcf.population import Population
from vcf.handler import Handler

ID_COLUMN = 2
INFO_COLUMN = 7
SUPERPOPULATIONS = {'AFR', 'AMR', 'EAS', 'EUR', 'SAS'}

__all__ = [
    'PredictorMariana',
    'PredictorJoao',
    'ALL_PREDICTORS',
]

class Predictor(Handler):
    """
    Bare bone implementation of a `Handler` that returns a prediciton for the
    population of individuals, based on the information in one VCF file
    """
    
    def __init__(self):
        super().__init__()
        
        self.population: Population
        self.polymorphisms: Set[str] = None
        self.polymorphism_count = 0
        
        # Initialize a dictionary that converts individual identifiers to
        # their index in the VCF file
        self.individuals: List[str]
        self.individual_indeces: Dict[str, int] = {}
        self.population_labels: List[str]
        self.from_info = False
        
        # Initialize an empty dictionary that will contain the monoid
        # structures for each individual
        self.structures: Dict[str, Any]
        
        self.labels: Dict[str, str] = {}
    
    
    def set_populations(self, population: Population, individuals: List[str]):
        self.population = population
        self.individuals = individuals
        
        self.from_info = set(self.population.groups()) == SUPERPOPULATIONS
    
    
    def set_polymorphisms(self, polymorphisms):
        self.polymorphisms = set(polymorphisms)
    
    
    def process_meta(self, line):
        pass
    
    
    def process_individuals(self, individuals):
        """
        Make sure all individuals given in the initialization method are present
        and create a monoid structure for each
        """
        
        if self.individuals is None:
            self.individuals = individuals
        
        self.population_labels = [
            self.population.individual_to_group[i]
            if self.population.has_individual(i) else '---'
            for i in individuals
        ]
        self.structures = {i: self.make_structure() for i in self.individuals}
        
        for idx, vcf_individual in enumerate(individuals):
            if vcf_individual in self.individuals:
                self.individual_indeces[vcf_individual] = idx
        
        if len(self.individuals) > len(self.individual_indeces):
            # Find the first that is not present
            for individual in self.individuals:
                if individual not in individuals:
                    raise Exception(
                        f'Individual {individual} not present in the VCF file'
                    )
    
    
    def process_variant(self, fixed_fields, format_field, variant_field):
        """
        Return exactly the GT data for a particular VCF field. For each
        individual, remove the non-mutated alleles, sort and remove duplicates.
        Equality is defined here as whether the resulting sets are equal.
        """
        
        #pylint: disable=W0613,R0201
        
        # Get the genotype of this individual
        genotype = variant_field.split(':')[0]
        
        if '|' in genotype:
            alleles = genotype.split('|')
        elif '/' in genotype:
            alleles = genotype.split('/')
        else:
            alleles = [genotype]
        
        alternatives = sum(1 for i in alleles if i != '0')
        return alternatives / len(alleles)
    
    
    def process_data(self, fixed_fields, format_field, data_fields):
        if self.polymorphisms is not None:
            if fixed_fields[ID_COLUMN] not in self.polymorphisms:
                return
            self.polymorphism_count += 1
        
        genotypes = [
            self.process_variant(fixed_fields, format_field, i)
            for i in data_fields
        ]
        
        freqs = self.get_frequencies(fixed_fields, genotypes)
        genotypes = self.filter_genotypes(genotypes)
        
        cache = {}
        
        for individual, genotype in zip(self.individuals, genotypes):
            if genotype in cache:
                value_list = cache[genotype]
            else:
                value_list = self.process_individual_genotype(freqs, individual, genotype, fixed_fields)
                cache[genotype] = value_list
            
            for value in value_list:
                self.update_individual(individual, value)
        
        if self.polymorphisms and self.polymorphism_count == len(self.polymorphisms):
            self.terminate_early = True
    
    
    def filter_genotypes(self, genotypes):
        """
        Return the genotypes of the relevant individuals.
        """
        
        return [genotypes[i] for i in self.individual_indeces.values()]
    
    
    def get_frequencies(self, fixed_fields, genotypes):
        """
        Either return the frequency for each population based on the alelles
        and population label for each individual, or extract the frequencies
        from the INFO column
        """
        
        if self.from_info:
            info = dict(
                i.split('=')
                for i in fixed_fields[INFO_COLUMN].split(';')
                if '=' in i)
            
            result = {}
            for pop in SUPERPOPULATIONS:
                if pop + '_AF' in info:
                    result[pop] = sum(float(i) for i in info[pop + '_AF'].split(','))
                else:
                    break
            
            if len(result) != len(SUPERPOPULATIONS):
                self.from_info = False
            else:
                return result
        
        counts = defaultdict(float)
        totals = defaultdict(int)
        for pop, f in zip(self.population_labels, genotypes):
            if pop == '---':
                continue
            
            counts[pop] += f
            totals[pop] += 1
        
        result = {pop: c/totals[pop] for pop, c in counts.items()}
        
        return result
    
    
    def process_individual_genotype(self, freqs, individual, genotype, fixed_fields):
        raise NotImplementedError
    
    
    def make_structure(self):
        """
        Create an empty monoid structure to hold on to the information of each
        individual.
        
        By default, the structure is a list.
        """
        
        #pylint: disable=R0201
        
        return []
    
    
    def append_to_structure(self, structure, value):
        """
        Append a value to a monoid structure.
        """
        
        #pylint: disable=R0201
        
        structure.append(value)
    
    
    def update_individual(self, individual, value):
        """
        Update the monoid structure of an individual with the provided value
        """
        
        self.append_to_structure(self.structures[individual], value)


class PredictorJoao(Predictor):
    
    def __init__(self, t1=0, t2=0):
        super().__init__()
        
        #pylint: disable=C0103
        
        self.t1 = t1
        self.t2 = t2
    
    def process_individual_genotype(self, freqs, individual, genotype, fixed_fields):
        d = {
            pop: abs(genotype - f)
            for pop, f in freqs.items()
        }
        ordered = sorted(d, key=lambda x: d[x])
        if d[ordered[1]] - d[ordered[0]] >= self.t1:
            return [ordered[0]]
        else:
            return []
    
    def terminate(self):
        relative = 0 < self.t2 < 1
        
        for individual, l in self.structures.items():
            ordered = Counter(l).most_common()
            diff = ordered[0][1] - ordered[1][1]
            
            if relative:
                diff /= ordered[0][1]
            
            if diff >= self.t2:
                self.labels[individual] = ordered[0][0]
            else:
                self.labels[individual] = '---'


class PredictorMariana(Predictor):
    
    def __init__(self, t1=0, t2=0):
        super().__init__()
        
        #pylint: disable=C0103
        
        self.t1 = t1
        self.t2 = t2
    
    def process_individual_genotype(self, freqs, individual, genotype, fixed_fields):
        result = []
        pairs = combinations(freqs, 2)
        for p1, p2 in pairs:
            d = abs(genotype - freqs[p1]) - abs(genotype - freqs[p2])
            if abs(d) >= self.t1:
                if d < 0:
                    result.append(p1)
                else:
                    result.append(p2)
        return result
    
    def terminate(self):
        relative = 0 < self.t2 < 1
        
        for individual, l in self.structures.items():
            ordered = Counter(l).most_common()
            diff = ordered[0][1] - ordered[1][1]
            
            if relative:
                diff /= ordered[0][1]
            
            if diff >= self.t2:
                self.labels[individual] = ordered[0][0]
            else:
                self.labels[individual] = '---'
        
ALL_PREDICTORS = {
    'PredictorJoao': PredictorJoao,
    'PredictorMariana': PredictorMariana,
}
