import csv
import math
import functools
from collections import defaultdict
from ortools.sat.python import cp_model
import sys
import argparse

try:
    import configparser
except ImportError:
    import ConfigParser


fiber_solver_defaults = {
        "threshold": 15,
        "core_penalty": 600,
        "max_cost": 4000,
        "max_chain": 2,
        "max_chain_core": 2
}


THRESHOLD = 15
CORE_PENALTY = 600
MAX_COST = 4000



class Connection(object):
    def __init__(self, cores,length):
        self.cores = int(cores)
        self.length = round(float(length))

    def is_multicore(self):
        return True if self.cores > 2 else False
    
    def __eq__(self,other):
        if issubclass(other,Connection):
            if self.cores == other.cores and self.length == other.length:
                return True
        return False
    
    def __len__(self):
        return 0 if self.cores <= 0 else self.length

    def __lt__(self, other):
        if issubclass(type(other),Connection):
            if self.cores < other.cores and self.length < other.length:
                return True
        return False
    
    def __le__(self, other):
        if issubclass(type(other),Connection):
            if self.cores <= other.cores and self.length < other.length:
                return True
        return False
    
    def __gt__(self, other):
        if issubclass(type(other),Connection):
            if self.cores > other.cores and self.length > other.length:
                return True
        return False
    
    def __gt__(self, other):
        if issubclass(type(other),Connection):
            if self.length >= other.length and self.cores >= other.cores:
                return True
        return False
    
    def __hash__(self):
        return hash((self.length, self.cores))


def cost(target,option):
    """
    Add penalty for if we overshoot the core count or we try to assign a huge multicore
    """
    penalty = 0
    if option.is_multicore() and not target.is_multicore() or target.cores != option.cores:
        penalty = 10000
    return penalty
    

class Fiber(Connection):
    def __init__(self, cores, length, name, armor):
        super(Fiber, self).__init__(cores, length)
        self.name = name
        self.armor = bool(armor)

    def __str__(self):
        str_core = "core"
        if self.cores > 1:
            str_core = str_core+"s"
        fib_info = ""
        if self.armor:
            fib_info = "Armored"
        return f"{self.name} ({self.length}m, {self.cores} {str_core} <{fib_info}>)"
   
    def __hash__(self):
          return hash(str(self))

class Link(Connection):
    def __init__(self, cores, length, source, dest):
        super(Link, self).__init__(cores, length)
        self.source = source
        self.dest = dest
    
    def __str__(self):
        str_core = "core"
        if self.cores > 1:
            str_core = str_core+"s"
        return f"{self.source} to {self.dest} ({self.length}m, {self.cores} {str_core})"
    
    def __hash__(self):
          return hash(str(self))

class Solver():
    @staticmethod
    def comp(fib1, fib2):
        if fib1["cores"] - fib2["cores"] == 0:
            return abs(fib1["length"] - fib2["length"])
        else:
            return abs(fib1["cores"] - fib2["cores"])

    @staticmethod
    def cost(target,option):
        """
        Add penalty for if we overshoot the core count or we try to assign a huge multicore
        """
        penalty = 0
        if target["cores"] != option["cores"] and (target["cores"] != 1 or option["cores"] > 3):
            penalty = 10000
        return penalty


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Assign fiber to links")
    parser.add_argument(
        "fibers",
        metavar="smf",
        type=argparse.FileType(),
        help="filename containing the fibers",
    )
    parser.add_argument(
        "links",
        metavar="links",
        type=argparse.FileType(),
        help="filename containing the links",
    )
    args = parser.parse_args()
    # Mark fiber with notes or broken cores/ need repairs as higher
    working_fibers = []
    reader = csv.DictReader(args.fibers, delimiter=",")
    for k in reader:
        working_fibers.append(Fiber(k["cores"], k["length"], k["name"],k["armor"]))

    needed_links = []
    reader = csv.DictReader(args.links, delimiter=",")
    for k in reader:
        if k["Subtype"].strip() == "Fibre":
            needed_links.append(Link(k["Cores"], k["Length"], k["From-Location"],k["To-Location"]))

    used_fibers = []

    model = cp_model.CpModel()

    x = {}

    # Populate dicitionary of assigments
    for fiber in working_fibers:
        for link in needed_links:
            x[fiber, link] = model.NewBoolVar(f"x[{fiber},{link}]")

    # Make sure we only assign each fiber once
    for fiber in working_fibers:
        model.AddAtMostOne(x[fiber, link] for link in needed_links)

    # Contraint fibers to only be assigned to links where they fit.
    # This requirement can be removed if running multiple fibres to reach
    # the amount of cores is acceptable
    for fiber in working_fibers:
        for link in needed_links:
            if fiber < link:
                model.Add(x[fiber, link] == False)

    # Optimize for the minize difference between fiber cores on chained fibers, or a new minimize target or a liniar contrainst
    for link in needed_links:
        # Make sure we do not daisy chain multi core fibers
        if link.is_multicore():
            model.Add(
                sum(x[fiber, link] for fiber in working_fibers) == 1
            )
            # Search Total length must be the length or 1.8 times more
            model.AddLinearConstraint(
                sum(len(fiber) * x[fiber, link] for fiber in working_fibers),
                len(link),
                int( len(link) * 1.8),
            )
        else:
            model.AddLinearConstraint(
                sum(len(fiber) * x[fiber, link] for fiber in working_fibers),
                len(link),
                len(link) * 5,
            )
            # We alloy a maximum of two fibers to be daisy chained
            model.AddLinearConstraint(
                sum(1 * x[fiber, link] for fiber in working_fibers),1,2)

    objective_terms = []

    costs = dict()

    for link in needed_links:
        for fiber in working_fibers:
            if not costs.get(link):
                costs[link] = dict()
            costs[link][fiber] = cost(link, fiber)

    for link in needed_links:
        for fiber in working_fibers:
            objective_terms.append(
                costs[link][fiber] * x[fiber, link]
            )

    # Minimise the amount of fibers used
    model.Minimize(
        sum(
            1000 * sum(x[fiber, link] for fiber in working_fibers)
            for link in needed_links
        )
    )

    for link in needed_links:
        # For each link minimize the total fiber count
        model.Minimize(
            sum(1000000 * x[fiber, link] for fiber in working_fibers)
        )
        # Minimize length overrun per link
        model.Minimize(
            sum(
                len(fiber) * x[fiber, link]
                for fiber in working_fibers
            )
            - len(link)
            + sum(1 * x[fiber, link] for fiber in working_fibers)
        )

        # Minimize core count overrun
        model.Minimize(
            sum(
                100 * fiber.cores * x[fiber, link]
                for fiber in working_fibers
            )
        )

    # Minimize the total core count over shoot penalty
    model.Minimize(sum(objective_terms))

    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    count_dict = defaultdict(int)
    coupler_count = defaultdict(int)
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        print(f"Total cost = {solver.ObjectiveValue()} \n")
        for link in needed_links:
            fib_list = []
            for fiber in working_fibers:
                if solver.BooleanValue(x[fiber, link]):
                    fib_list.append(fiber)
                    # print(f'Fiber {fiber["name"]} assigned to link {link["name"]}.' +
                    # f' {fiber} {link}')
            fib_names = ", ".join(list(map(lambda x: str(x), fib_list)))
            core_list = (list(map(lambda x: x.cores, fib_list)))
            max_core = max(core_list)
            
            couplersize = 1
            
            if max_core >= 4:
                couplersize = 4
            elif max_core >= 2:
                couplersize = 2
            
            if len(core_list) > 1:
                coupler_count[couplersize] = coupler_count[couplersize] + (len(core_list)-1 * (max_core/couplersize)) 
            print(
                f'Link {link} has [ {fib_names} ]'
            )
            count_dict[max_core] = count_dict[max_core] + 1
    else:
        print("No solution found.")
    for fiber in working_fibers:
        used = False
        for link in needed_links:
            if solver.BooleanValue(x[fiber, link]):
                used = True
        if used == False:
            print(f"not using {fiber}")

    # TODO make a matrix of all allowed combinations and match on that

    for key, value in count_dict.items():
        print(f"cores {key} amount {value}")

    min_coupler = 2 
    for key, value in coupler_count.items():
        if key < min_coupler:
            coupler_count[key] = coupler_count[min_coupler]+ value
    
    for key, value in coupler_count.items():
        if key < min_coupler:
            continue
        print(f"coupler {key} amount {value}")
