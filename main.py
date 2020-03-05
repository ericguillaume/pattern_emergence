import random

import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import dgamma, norm
from tqdm import tqdm


# todo(eric) add food instability
# todo(eric) reproduce with random gene
# groups seems to exterminate each  other too well already
# @ punir un group trop present par virus !!! plus de chance pour eux de le choper
# merge les cells en une seul (u on peut avoir jusqu a 10 animauxc etc..., plus simple)


'''
    new born should have same energy, maybe there is more
    
    reproduction with only 0.5 given to new baby increases diversity destruction,
        maximiser energy n'y change rien, energie donner par food bien gausianniser n'y change rien  !!
        
    todo(eric) coder 10 animaux... ont tous chances de crever.. puis un est copier ave cchance... voir ce que ca donne
    
    todo(eric) coder premiers examples 2/3 iterations
        
        
        
        
    SCENARIOS: size 20      

    no viruses, IDENTICAL reproduction, 20 animals start, animals dont walk on each others 
        => diversity kept (for 800 steps periods only)
        
    no viruses, IDENTICAL reproduction, 20 animals start, 10k steps
        => diversity destroyed 
        
    no viruses, IDENTICAL reproduction, 20 animals start, no more SOLIDARITY, 10k steps
        => diversity destroyed 
        
    no viruses, IDENTICAL reproduction, 20 animals start, no more SOLIDARITY, no more geographical BOUNDARIES,
            10k steps
        => diversity destroyed 
    no viruses, IDENTICAL reproduction, 20 animals start, no more SOLIDARITY, no more geographical BOUNDARIES,
            10k steps, reproduced next, no move move
        => diversity destroyed , 400 animals, energy per animals 1.1 (same as before)
        
    10 familles de 200 membres => extinction ralentie
'''

random.seed(8)


class Gene:
    def __init__(self, value):
        self.value = value

    @staticmethod
    def new_random_gene():
        value = random.uniform(0.0, 1.0)
        return Gene(value)

    def reproduce_identical(self):
        return Gene(self.value)


class Animal:
    DEFAULT_ENERGY = 1.0
    REPRODUCTION_ENERGY = DEFAULT_ENERGY / 3.0
    MAX_ENERGY = DEFAULT_ENERGY * 2.0

    def __init__(self, x, y, energy, is_in_control_group, gene):
        self.is_in_control_group = is_in_control_group
        self.gene = gene
        self.x = x
        self.y = y
        self.energy = energy

    @staticmethod
    def new_default_energy_animal(x, y):
        gene_value = int(random.uniform(0.0, 10.0)) / 10.0
        gene = Gene(gene_value)

        # is_in_control_group = (random.randint(0, 1) == 0)
        is_in_control_group = gene_value < 0.5
        # gene_value = 0.0 if is_in_control_group else 1.0
        # gene = Gene(gene_value)
        # gene = Gene.new_random_gene()
        return Animal(x, y, Animal.DEFAULT_ENERGY, is_in_control_group, gene)

    def move(self, universe_size, move_size):
        new_x = (self.x + random.randint(-move_size, move_size)) % universe_size
        new_y = (self.y + random.randint(-move_size, move_size)) % universe_size

        self.x = new_x
        self.y = new_y
        return self

    def get_coordinates(self):
        return self.x, self.y

    def add_energy(self, amount):
        self.energy = min(Animal.MAX_ENERGY, self.energy + amount)

    def consume_energy(self, amount):
        self.energy -= amount

    def is_alive(self):
        return self.energy >= 0.

    def reproduce(self, universe_size):
        '''
        :return: new animal
        '''
        self.consume_energy(Animal.REPRODUCTION_ENERGY)
        return Animal(self.x, self.y, Animal.REPRODUCTION_ENERGY, self.is_in_control_group,
                      self.gene.reproduce_identical()) \
            .move(universe_size, 3)


class Universe:
    ANIMALS = 2000  # 6
    SIZE = 80  # 20  # 12
    MOVE_SIZE = 2  # 2
    REPRODUCTION_MIN_ENERGY = 1.5
    ENERGY_CONSUMED_PER_STEP = 0.05
    ENERGY_LOST_BY_WALKING_ONTO_EACH_OTHER = ENERGY_CONSUMED_PER_STEP
    POSITIVE_FOOD_INPUT = ENERGY_CONSUMED_PER_STEP / 5

    def __init__(self):
        self.animals = set([self.new_random_animal() for _ in range(Universe.ANIMALS)])

    def run_step(self):
        # help each others (control_group only)
        # self.animals_help_each_others(filter_out_control_group=True)

        # walk on each others
        self.animals_walk_on_each_others()

        # energy loss and death
        deads = []
        for ani in self.animals:
            ani.consume_energy(Universe.ENERGY_CONSUMED_PER_STEP)
            if not ani.is_alive():
                deads.append(ani)
        for dead in deads:
            self.animals.remove(dead)

        # food
        for ani in self.animals:
            x, y = ani.get_coordinates()
            food = self.get_food(x, y)
            ani.add_energy(food)

        # move
        for ani in self.animals:
            ani.move(Universe.SIZE, Universe.MOVE_SIZE)

        # reproduction
        new_animals = []
        for ani in self.animals:
            if ani.energy > Universe.REPRODUCTION_MIN_ENERGY:
                new_animals.append(ani.reproduce(Universe.SIZE))
        for new_ani in new_animals:
            self.animals.add(new_ani)

        energy_control = sum([ani.energy for ani in self.animals if ani.is_in_control_group])
        energy_test = sum([ani.energy for ani in self.animals if not ani.is_in_control_group])
        animals_control = len([ani for ani in self.animals if ani.is_in_control_group])
        animals_test = len([ani for ani in self.animals if not ani.is_in_control_group])
        genes = [ani.gene.value for ani in self.animals]
        return {"animals_control": animals_control,
                "animals_test": animals_test,
                "energy_control": energy_control,
                "energy_test": energy_test}, genes

    def animals_walk_on_each_others(self):
        locations_animals_dict = self.build_locations_seq_animals_dict()
        for _, animals in locations_animals_dict.items():
            if len(animals) == 1:
                continue
            for ani in animals:
                ani.consume_energy(Universe.ENERGY_LOST_BY_WALKING_ONTO_EACH_OTHER)

    def build_locations_seq_animals_dict(self):
        locations_animals_dict = {}
        for ani in self.animals:
            x_y = ani.get_coordinates()
            if x_y not in locations_animals_dict:
                locations_animals_dict[x_y] = []
            locations_animals_dict[x_y].append(ani)
        return locations_animals_dict

    def animals_help_each_others(self, filter_out_control_group):
        '''
        it doesn't manage 2 animals on the same location
        '''
        locations_animals_dict = self.build_locations_animals_dict(
            filter_out_control_group=filter_out_control_group)
        for (x, y), ani_receiver in locations_animals_dict.items():
            if ani_receiver.energy > 0.1:
                continue
            for other_x in range(x - Universe.SOLIDARITY_LENGTH, x + Universe.SOLIDARITY_LENGTH + 1):
                if other_x < 0 or other_x >= Universe.SIZE:
                    continue
                for other_y in range(y - Universe.SOLIDARITY_LENGTH, y + Universe.SOLIDARITY_LENGTH + 1):
                    if other_y < 0 or other_y >= Universe.SIZE:
                        continue
                    if not (other_x, other_y) in locations_animals_dict:
                        continue
                    ani_giver = locations_animals_dict[(other_x, other_y)]
                    if ani_giver.energy > Universe.SOLIDARITY_MIN_ENERGY_TO_HELP:
                        ani_giver.consume_energy(Universe.SOLIDARITY_GIFT + Universe.SOLIDARITY_COST)
                        ani_receiver.add_energy(Universe.SOLIDARITY_GIFT)

    def build_locations_animals_dict(self, filter_out_control_group=False):
        locations_animals_dict = {}
        for ani in self.animals:
            if filter_out_control_group and ani.is_in_control_group:
                continue
            x, y = ani.get_coordinates()
            locations_animals_dict[(x, y)] = ani
        return locations_animals_dict

    def get_food(self, x, y):
        loc = Universe.POSITIVE_FOOD_INPUT + Universe.ENERGY_CONSUMED_PER_STEP
        scale = Universe.ENERGY_CONSUMED_PER_STEP / 10.0
        result = norm.rvs(loc=loc, scale=scale, size=1)[0]
            # random.uniform(Universe.POSITIVE_FOOD_INPUT,
            #                   Universe.POSITIVE_FOOD_INPUT + 2 * Universe.ENERGY_CONSUMED_PER_STEP)
        result = loc
        return max(result, 0.0)

    @staticmethod
    def new_random_animal():
        x, y = Universe.new_random_coordinates()
        return Animal.new_default_energy_animal(x, y)

    @staticmethod
    def new_random_coordinates():
        x = random.randint(0, Universe.SIZE - 1)
        y = random.randint(0, Universe.SIZE - 1)
        return x, y


def display_hist_genes_values():
    plt.hist(genes_values, bins=100)
    plt.title("Gene values")
    plt.show()


STEPS = 10000
GENES_HIST_DISPLAYED = 0  # 10
GENES_HIST_DISPLAYED_AT_END = True

uni = Universe()
analytics = []
for idx, _ in enumerate(tqdm(range(STEPS))):
    row_analytics, genes_values = uni.run_step()
    analytics.append(row_analytics)
    if GENES_HIST_DISPLAYED > 0 and idx % int(STEPS / GENES_HIST_DISPLAYED) == 0:
        display_hist_genes_values()

if GENES_HIST_DISPLAYED_AT_END:
    display_hist_genes_values()


analytics = pd.DataFrame(analytics, columns=["animals_control", "animals_test", "energy_control", "energy_test"])
analytics["animals"] = analytics["animals_control"] + analytics["animals_test"]
analytics[["animals_control", "animals_test", "animals"]].plot()
plt.show()

analytics["energy"] = analytics["energy_control"] + analytics["energy_test"]
analytics["energy_per_animal"] = analytics["energy"] / analytics["animals"]
analytics[["energy_per_animal"]].plot()
plt.show()
