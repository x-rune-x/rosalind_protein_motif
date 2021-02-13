import requests


# Create a class that links each protein ID to its specific sequence and allow for operations on the sequence.
class ProteinFASTA:
    def __init__(self, fasta_id, sequence):
        self.fasta_id = fasta_id
        self.sequence = sequence

    def get_id(self):
        return self.fasta_id

    def get_seq(self):
        return self.sequence

    def _get_length(self):
        return len(self.sequence)

    # Finds the starting position of any N-glycosylation motifs in the sequence and stores their starting position.
    def find_motif(self):
        seq_length = len(self.sequence)
        pos_list = []

        for position in range(seq_length):
            frame = self.sequence[position: position + 4]
            if len(frame) == 4:
                # N-glycosylation motif shorthand is N{P}[ST]{P}
                if frame[0] == "N" and frame[1] != "P" and (frame[2] == "S" or frame[2] == "T") and frame[3] != "P":
                    # Add 1 to position because Python index begins at 0 while standard protein length begins at 1.
                    pos_list.append(position + 1)

        return pos_list


# Take a list of protein IDs and retrieve their specific sequence from UniProt.
def get_seqs(id_file):
    id_list = open(id_file, "r")
    fasta_list = []

    for line in id_list:
        protein_id = line.rstrip()
        url = f"http://www.uniprot.org/uniprot/{protein_id}.fasta"
        sequence = requests.get(url)
        id_txt = sequence.text
        line_list = id_txt.splitlines()
        seq_txt = ""
        for seq_line in range(1, len(line_list)):
            seq_txt += line_list[seq_line]

        fasta_list.append(ProteinFASTA(protein_id, seq_txt))

    id_list.close()
    return fasta_list


def find_positions(input_list):
    out_file = open("positions.txt", "w")
    for line in input_list:
        if len(line.find_motif()) != 0:
            out_file.write(line.get_id() + "\n")
            positions = ""
            for position in line.find_motif():
                positions += str(position) + " "
            out_file.write(positions + "\n")

    out_file.close()


sequence_list = get_seqs("rosalind_mprt.txt")
find_positions(sequence_list)
