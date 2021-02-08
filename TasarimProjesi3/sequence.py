class Sequence:
    @staticmethod
    def read(filename):
        names = []
        sequences = []
        with open(filename) as sequence_file:
            lines = sequence_file.readlines()

            for i in range(0, len(lines) - 1, 2):
                names.append(lines[i].strip())
                sequences.append(lines[i + 1].strip())

            return names, sequences

