import os

if __name__ == "__main__":
    with open(os.path.join(os.path.dirname(__file__), "recipe.yaml")) as f:
        lines = f.readlines()
        for i in range(len(lines) - 1):
            if lines[i].strip() == "build:" and lines[i + 1].strip().startswith("number:"):
                lines[i + 1] = f"  number: {int(lines[i + 1].strip().split()[-1]) + 1}\n"
                break

    with open(os.path.join(os.path.dirname(__file__), "recipe.yaml"), "w") as f:
        f.writelines(lines)
