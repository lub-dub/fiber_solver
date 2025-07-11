import csv
import math
import functools
from collections import defaultdict
from ortools.sat.python import cp_model
import sys
import argparse

import tomllib


class Connection(object):
    def __init__(self, cores, length):
        self.cores = int(cores)
        self.length = round(float(length))

    def is_multicore(self):
        return True if self.cores > 2 else False

    def __eq__(self, other):
        if issubclass(other, Connection):
            if self.cores == other.cores and self.length == other.length:
                return True
        return False

    def __len__(self):
        return 0 if self.cores <= 0 else self.length

    def __lt__(self, other):
        if issubclass(type(other), Connection):
            if self.cores < other.cores and self.length < other.length:
                return True
        return False

    def __le__(self, other):
        if issubclass(type(other), Connection):
            if self.cores <= other.cores and self.length < other.length:
                return True
        return False

    def __gt__(self, other):
        if issubclass(type(other), Connection):
            if self.cores > other.cores and self.length > other.length:
                return True
        return False

    def __gt__(self, other):
        if issubclass(type(other), Connection):
            if self.length >= other.length and self.cores >= other.cores:
                return True
        return False

    def __hash__(self):
        return hash((self.length, self.cores))


def cost(target, option):
    """
    Add penalty for if we overshoot the core count or we try to assign a huge multicore
    """
    penalty = 0
    if (
        option.is_multicore()
        and not target.is_multicore()
        or target.cores != option.cores
    ):
        penalty = 10000
    return penalty


class Fiber(Connection):
    def __init__(self, cores, length, name, armor):
        super(Fiber, self).__init__(cores, length)
        self.name = name
        self.armor = bool(int(armor))

    def __str__(self):
        str_core = "core"
        fib_info = []
        if self.cores > 1:
            str_core = str_core + "s"
        if self.cores > 2:
            fib_info.append("Multicore")
        if self.armor:
            fib_info.append("Armored")
        return f"{self.name} ({self.length}m, {self.cores} {str_core} <{','.join(fib_info)}>)"

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
            str_core = str_core + "s"
        return f"{self.source} to {self.dest} ({self.length}m, {self.cores} {str_core})"

    def __hash__(self):
        return hash(str(self))


class Solver:

    def __init__(self, config=None):
        default_config = {
            "slack": 5,
            "core_penalty": 10000,
            "max_cost": 4000,
            "max_length_multicore": 1.8,
            "max_length": 5,
            "max_chain": 2,
            "max_cores":24,
            "max_chain_core": 2,
            "min_coupler": 2,
            "max_coupler": 4,
        }
        self.config = config if config else default_config

    def load(self, links, fibers):
        self.model = cp_model.CpModel()
        self.x = {}
        self.costs = {}
        self.objective_terms = []
        self.status = None
        self.solver = None
        self.fibers = fibers
        self.links = links

        for fiber in fibers:
            for link in links:
                self.x[fiber, link] = self.model.NewBoolVar(f"x[{fiber},{link}]")
                if fiber < link:
                    self.model.Add(self.x[fiber, link] == False)
        # every fiber only once
        for fiber in fibers:
            self.model.AddAtMostOne(self.x[fiber, link] for link in links)

        for link in links:
            # Make sure we do not daisy chain multi core fibers
            if link.is_multicore():
                self.model.Add(sum(self.x[fiber, link] for fiber in fibers) == 1)
                # Search Total length must be the length or 1.8 times more
                self.model.AddLinearConstraint(
                    sum(len(fiber) * self.x[fiber, link] for fiber in fibers),
                    len(link),
                    int(len(link) * self.config["max_length_multicore"]),
                )
                self.model.AddLinearConstraint(
                    sum(fiber.cores * self.x[fiber, link] for fiber in fibers),
                    link.cores,
                    int(self.config["max_cores"]),
                )
            else:
                self.model.AddLinearConstraint(
                    sum(len(fiber) * self.x[fiber, link] for fiber in fibers) + self.config["slack"],
                    len(link),
                    len(link) * self.config["max_length"],
                )
                # We alloy a maximum of two fibers to be daisy chained
                self.model.AddLinearConstraint(
                    sum(1 * self.x[fiber, link] for fiber in working_fibers), 1, self.config["max_chain"]
                )
        for link in links:
            for fiber in fibers:
                if not self.costs.get(link):
                    self.costs[link] = dict()
                self.costs[link][fiber] = self.cost(link, fiber)
                self.objective_terms.append(self.costs[link][fiber] * self.x[fiber, link])

        self.model.Minimize(
            sum(1000 * sum(self.x[fiber, link] for fiber in fibers) for link in links)
        )

        for link in links:
            # For each link minimize the total fiber count
            self.model.Minimize(sum(50000 * self.x[fiber, link] for fiber in fibers))
            # Minimize length overrun per link
            self.model.Minimize(
                sum(len(fiber) * self.x[fiber, link] for fiber in fibers) - len(link)
            )

            # Minimize core count overrun
            self.model.Minimize(
                sum(10000 * fiber.cores * self.x[fiber, link] for fiber in fibers)
            )

        # Minimize the total core count over shoot penalty
        self.model.Minimize(sum(self.objective_terms))

    def solve(self):
        self.solver = cp_model.CpSolver()
        self.status = self.solver.Solve(self.model)

    def __str__(self):
        if self.solver == None:
            return "Not yet solved"

        if self.status != cp_model.OPTIMAL and status != cp_model.FEASIBLE:
            return "No solution"

        used_fibers = []
        count_dict = defaultdict(int)
        coupler_count = defaultdict(int)
        ret = ""
        ret += f"Total cost = {self.solver.ObjectiveValue()} \n"
        for link in self.links:
            fib_list = []
            for fiber in self.fibers:
                if self.solver.BooleanValue(self.x[fiber, link]):
                    fib_list.append(fiber)
            fib_names = ", ".join(list(map(lambda x: str(x), fib_list)))
            core_list = list(map(lambda x: x.cores, fib_list))
            max_core = max(core_list)

            couplersize = 1

            if max_core >= self.config["max_coupler"]:
                couplersize = self.config["max_coupler"]
            elif max_core >= self.config["min_coupler"]:
                couplersize = self.config["min_coupler"]

            if len(core_list) > 1:
                coupler_count[couplersize] = coupler_count[couplersize] + (
                    len(core_list) - 1 * (max_core / couplersize)
                )
            ret +=f"Link {link} has [ {fib_names} ]\n"
            count_dict[max_core] = count_dict[max_core] + 1
        for fiber in self.fibers:
            used = False
            for link in self.links:
                if self.solver.BooleanValue(self.x[fiber, link]):
                    used = True
            if used == False:
                ret +=f"not using {fiber}\n"

        for key, value in sorted(count_dict.items(), key=lambda x: x[0]):
            ret += f"cores {key} amount {value}\n"

        # minimum coupler size, so if we couple simplex cable but only have duplex couples increase the count
        min_coupler = self.config["min_coupler"]
        for key, value in coupler_count.items():
            if key < min_coupler:
                coupler_count[key] = coupler_count[min_coupler] + value

        for key, value in coupler_count.items():
            if key < min_coupler:
                continue
            ret += f"coupler {key} amount {value}\n"
        return ret

    def cost(self, target, option):
        """
        Add penalty for if we overshoot the core count or we try to assign a huge multicore
        """
        penalty = 0
        if (
            option.is_multicore()
            and not target.is_multicore()
            or target.cores != option.cores
        ):
            penalty = self.config["core_penalty"]
        return penalty


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Assign fiber to links")
    parser.add_argument(
        "--config", type=argparse.FileType("rb"), help="Filename with config"
    )
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
    if args.config:
        config = tomllib.load(args.config)
    # Mark fiber with notes or broken cores/ need repairs as higher
    working_fibers = []
    reader = csv.DictReader(args.fibers, delimiter=",")
    for k in reader:
        working_fibers.append(Fiber(k["cores"], k["length"], k["name"], k["armor"]))

    needed_links = []
    reader = csv.DictReader(args.links, delimiter=",")
    for k in reader:
        if k["Subtype"].strip() == "Fibre":
            needed_links.append(
                Link(k["Cores"], k["Length"], k["From-Location"], k["To-Location"])
            )
    solver = Solver()
    solver.load(needed_links,working_fibers)
    solver.solve()
    print(solver)
