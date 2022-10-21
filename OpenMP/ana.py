import hatchet as ht

if __name__ == "__main__":
    dirname = "hpctoolkit-sim16-database"
    gf = ht.GraphFrame.from_hpctoolkit(dirname)

print(gf.dataframe)

print(gf.tree(metric_column="time (inc)"))

with open("test.dot", "w") as dot_file:
    dot_file.write(gf.to_dot())
