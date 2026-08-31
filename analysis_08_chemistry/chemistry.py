import csv
from pathlib import Path

FAMILIES = ("H2O", "OH", "H2", "H", "CO2", "CO", "O2", "O", "C")


def reactions(path, f107p):
    with Path(path).open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 35 or f107p <= 0:
        raise ValueError("expected 35 reactions and positive F10.7P")
    for row in rows:
        row["A0"], row["A1"] = float(row["A0"]), float(row["A1"])
        row["rate_s-1"] = row["A0"] * f107p ** row["A1"]
        for name in FAMILIES:
            row[name] = float(row[name])
    return rows


def network(rows, cross_sections):
    result = {}
    for parent in FAMILIES:
        selected = [row for row in rows if row["parent"] == parent]
        total = sum(row["rate_s-1"] for row in selected)
        if total <= 0:
            raise ValueError(f"missing destruction rate for {parent}")
        children = []
        for child in FAMILIES:
            value = sum(row["rate_s-1"] * row[child] for row in selected) / total
            if value > 0:
                children.append((child, value))
        sigma = cross_sections.get(parent, cross_sections["daughter_proxy"])
        result[parent] = {
            "rate_1au_s-1": total,
            "sigma_cm2": sigma,
            "children": tuple(children),
        }
    return result


def velocity_key(parent, child):
    if parent == "H2O":
        return {"OH": "OH", "O": "O_H2O", "H": "H", "H2": "H2_photo"}[child]
    if parent == "OH":
        return {"O": "O_OH", "H": "H"}[child]
    if parent == "H2":
        return "H"
    return "heavy"


def chains(primary, network):
    output = []

    def visit(family, speed, segments, weight):
        current = segments + ((family, speed),)
        output.append((current, weight))
        if len(current) == 3:
            return
        for child, branch in network[family]["children"]:
            visit(child, velocity_key(family, child), current, weight * branch)

    visit(primary, primary, (), 1.0)
    return output
