import csv
import math
import functools
from collections import defaultdict
from ortools.sat.python import cp_model
import sys
import argparse


THRESHOLD = 15
CORE_PENALTY = 600
MAX_COST = 4000


def comp(fib1, fib2):
    if fib1["cores"] - fib2["cores"] == 0:
        return abs(fib1["length"] - fib2["length"])
    else:
        return abs(fib1["cores"] - fib2["cores"])


def cost(fib1, fib2):
    """
    Add penalty for if we overshoot the core count or we try to assign a huge multicore
    """
    penalty = 0
    if fib1["cores"] != fib2["cores"] and (fib1["cores"] != 1 or fib2["cores"] > 3):
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
        working_fibers.append(
            {
                "cores": int(k["cores"]),
                "armor": k["armor"],
                "length": int(k["length"]),
                "name": k["name"],
            }
        )

    needed_links = []
    reader = csv.DictReader(args.links, delimiter=",")
    for k in reader:
        if k["Subtype"].strip() == "Fibre":
            needed_links.append(
                {
                    "cores": int(k["Cores"]),
                    "length": math.ceil(float(k["Length"]) + 20),
                    "name": "{} to {}".format(k["From-Location"], k["To-Location"]),
                }
            )

    # Sort runs from most cores to least cores ie higher contraint fibers first
    needed_links = sorted(
        needed_links, key=lambda x: x["cores"] * 1000 + x["length"], reverse=True
    )

    used_fibers = []

    model = cp_model.CpModel()

    x = {}

    # Populate dicitionary of assigments
    for fiber in working_fibers:
        for link in needed_links:
            x[fiber["name"], link["name"]] = model.NewBoolVar(f"x[{fiber},{link}]")

    # Make sure we only assign each fiber once
    for fiber in working_fibers:
        model.AddAtMostOne(x[fiber["name"], link["name"]] for link in needed_links)

    # Contraint fibers to only be assigned to links where they fit.
    # This requirement can be removed if running multiple fibres to reach
    # the amount of cores is acceptable
    for fiber in working_fibers:
        for link in needed_links:
            if fiber["cores"] < link["cores"]:
                model.Add(x[fiber["name"], link["name"]] == False)

    # Optimize for the minize difference between fiber cores on chained fibers, or a new minimize target or a liniar contrainst
    for link in needed_links:
        # Make sure we do not daisy chain multi core fibers
        if link["cores"] > 2:
            model.Add(
                sum(x[fiber["name"], link["name"]] for fiber in working_fibers) == 1
            )
            # Search Total length must be the length or 1.8 times more
            model.AddLinearConstraint(
                sum(
                    fiber["length"] * x[fiber["name"], link["name"]]
                    for fiber in working_fibers
                ),
                link["length"],
                int(link["length"] * 1.8),
            )
        else:
            model.AddLinearConstraint(
                sum(
                    fiber["length"] * x[fiber["name"], link["name"]]
                    for fiber in working_fibers
                ),
                link["length"],
                int(link["length"] * 5),
            )
            # We alloy a maximum of two fibers to be daisy chained
            model.AddLinearConstraint(
                sum(1 * x[fiber["name"], link["name"]] for fiber in working_fibers),
                1,
                2,
            )

    objective_terms = []

    costs = dict()

    for link in needed_links:
        for fiber in working_fibers:
            if not costs.get(link["name"]):
                costs[link["name"]] = dict()
            costs[link["name"]][fiber["name"]] = cost(link, fiber)

    for link in needed_links:
        for fiber in working_fibers:
            objective_terms.append(
                costs[link["name"]][fiber["name"]] * x[fiber["name"], link["name"]]
            )

    # Minimise the amount of fibers used
    model.Minimize(
        sum(
            1000 * sum(x[fiber["name"], link["name"]] for fiber in working_fibers)
            for link in needed_links
        )
    )

    for link in needed_links:
        # For each link minimize the total fiber count
        model.Minimize(
            sum(1000000 * x[fiber["name"], link["name"]] for fiber in working_fibers)
        )
        # Minimize length overrun per link
        model.Minimize(
            sum(
                fiber["length"] * x[fiber["name"], link["name"]]
                for fiber in working_fibers
            )
            - link["length"]
            + sum(1 * x[fiber["name"], link["name"]] for fiber in working_fibers)
        )
        # Minimize core count overrun
        model.Minimize(
            sum(
                100 * fiber["cores"] * x[fiber["name"], link["name"]]
                for fiber in working_fibers
            )
        )

    # Minimize the total core count over shoot penalty
    model.Minimize(sum(objective_terms))

    solver = cp_model.CpSolver()
    status = solver.Solve(model)

    count_dict = defaultdict(int)

    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        print(f"Total cost = {solver.ObjectiveValue()} \n")
        for link in needed_links:
            fib_list = []
            for fiber in working_fibers:
                if solver.BooleanValue(x[fiber["name"], link["name"]]):
                    fib_list.append(fiber)
                    # print(f'Fiber {fiber["name"]} assigned to link {link["name"]}.' +
                    # f' {fiber} {link}')
            fib_names = list(map(lambda x: x["name"], fib_list))
            fib_length = list(map(lambda x: x["length"], fib_list))
            fib_cores = list(map(lambda x: x["cores"], fib_list))
            print(
                f'Link {link["name"]} with length {link["length"]} and {link["cores"]} cores has {fib_names} {fib_length} {sum(fib_length)} {fib_cores}'
            )
            count_dict[max(fib_cores)] = count_dict[max(fib_cores)] + 1
    else:
        print("No solution found.")
    for fiber in working_fibers:
        used = False
        for link in needed_links:
            if solver.BooleanValue(x[fiber["name"], link["name"]]):
                used = True
        if used == False:
            print(f"not using {fiber}")

    # TODO make a matrix of all allowed combinations and match on that

    for key, value in count_dict.items():
        print(f"cores {key} amount {value}")
