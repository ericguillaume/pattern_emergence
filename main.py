import random
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

random.seed(8)


class Gene:

    def __init__(self):
        self.value = random.uniform(0.0, 1.0)


class Animal:
    DEFAULT_ENERGY = 1.0

    def __init__(self, x, y, energy):
        self.is_in_control_group = (random.randint(0, 1) == 0)
        self.gene = Gene()
        self.x = x
        self.y = y
        self.energy = energy

    @staticmethod
    def new_default_animal(x, y):
        return Animal(x, y, Animal.DEFAULT_ENERGY)

    def move(self, universe_size, move_size):
        new_x = self.x + random.randint(-move_size, move_size)
        new_y = self.y + random.randint(-move_size, move_size)

        new_x = min(max(0, new_x), universe_size - 1)
        new_y = min(max(0, new_y), universe_size - 1)

        self.x = new_x
        self.y = new_y
        return self

    def get_coordinates(self):
        return self.x, self.y

    def add_energy(self, amount):
        self.energy += amount

    def consume_energy(self, amount):
        self.energy -= amount

    def is_alive(self):
        return self.energy >= 0.

    def reproduce(self, universe_size):
        '''
        :return: new animal
        '''
        self.energy /= 2
        return Animal(self.x, self.y, self.energy).move(universe_size, 1)


class Universe:
    ANIMALS = 6
    SIZE = 12
    MOVE_SIZE = 2
    REPRODUCTION_MIN_ENERGY = 2.0
    ENERGY_CONSUMED_PER_STEP = 0.05
    ENERGY_LOST_BY_WALKING_ONTO_EACH_OTHER = ENERGY_CONSUMED_PER_STEP
    POSITIVE_FOOD_INPUT = ENERGY_CONSUMED_PER_STEP / 5
    SOLIDARITY_MIN_ENERGY_TO_HELP = 1.0
    SOLIDARITY_LENGTH = 4
    SOLIDARITY_GIFT = 0.1
    SOLIDARITY_COST = ENERGY_CONSUMED_PER_STEP / 2

    def __init__(self):
        self.animals = set([self.new_random_animal() for _ in range(Universe.ANIMALS)])

    def run_step(self):
        # help each others
        self.animals_help_each_others()

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
        return {"animals_control": animals_control,
                "animals_test": animals_test,
                "energy_control": energy_control,
                "energy_test": energy_test}

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

    def animals_help_each_others(self):
        '''
        it doesn't manage 2 animals on the same location
        '''
        locations_animals_dict = self.build_locations_animals_dict(filter_is_in_group=True)
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

    def build_locations_animals_dict(self, filter_is_in_group=False):
        locations_animals_dict = {}
        for ani in self.animals:
            if filter_is_in_group and ani.is_in_control_group:
                continue
            x, y = ani.get_coordinates()
            locations_animals_dict[(x, y)] = ani
        return locations_animals_dict

    def get_food(self, x, y):
        return random.uniform(Universe.POSITIVE_FOOD_INPUT,
                              Universe.POSITIVE_FOOD_INPUT + 2 * Universe.ENERGY_CONSUMED_PER_STEP)

    @staticmethod
    def new_random_animal():
        x, y = Universe.new_random_coordinates()
        return Animal.new_default_animal(x, y)

    @staticmethod
    def new_random_coordinates():
        x = random.randint(0, Universe.SIZE - 1)
        y = random.randint(0, Universe.SIZE - 1)
        return x, y


uni = Universe()
analytics = []
for _ in tqdm(range(10000)):
    row_analytics = uni.run_step()
    analytics.append(row_analytics)

# analytics = pd.DataFrame(analytics, columns=["animals", "energy"])
# analytics["energy_per_animal"] = analytics["energy"] / analytics["animals"]
#
# analytics[["animals", "energy"]].plot()
# plt.show()
# analytics[["energy_per_animal"]].plot()
# plt.show()

analytics = pd.DataFrame(analytics, columns=["animals_control", "animals_test", "energy_control", "energy_test"])
analytics["animals"] = analytics["animals_control"] + analytics["animals_test"]
analytics[["animals_control", "animals_test", "animals"]].plot()
plt.show()
