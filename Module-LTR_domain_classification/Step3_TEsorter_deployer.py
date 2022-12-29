import os
import re
import pandas as pd


# Env req: /data5/condaenvs/AndroTE

class IndividualAccession(object):
    def __init__(self):
        # Root folder for the project, containing all links.
        self._accession_name = os.getcwd().split("/")[-1]  # GCA000000.1
        self._current_dir = os.getcwd() + "/"

    def split_TElib(self):
        frag_dict = {}  # TE_tag => TE_raw_sequence
        LTR_RT_comp = []  # [[TE1_tag, TE1_raw_sequence], [TE2_tag, TE2_raw_sequence]]
        non_LTR_RT_comp = []  # [[TE1_tag, TE1_raw_sequence], [TE2_tag, TE2_raw_sequence]]

        # Catalog FA entries
        current_frag_list = []
        frag_tag = ""
        with open(self._current_dir + self._accession_name + "_telib.fa", mode='r') as lines:
            for line in lines:
                if line.startswith(">"):
                    if current_frag_list:  # Handle The Fragment Before
                        frag_dict[frag_tag] = "".join(current_frag_list)
                        current_frag_list.clear()

                    frag_tag = line.rstrip("\n")
                    continue

                current_frag_list.append(line.rstrip("\n"))
            frag_dict[frag_tag] = "".join(current_frag_list)  # Handle The Last Fragment
        print("Parse EDTA TElib completed, identified %d entries." % (len(frag_dict.keys())))

        # Parse and filter LTR-RT entries
        match_rule_LTR_RT = re.compile(r'^>(?P<TEID>TE_\d+)_(LTR|INT)#LTR/(?P<SUPFAM>Gypsy|Copia|unknown)$')
        copia_count = 0
        gypsy_count = 0
        unknown_count = 0
        for k, v in frag_dict.items():
            TE_clade = re.search(match_rule_LTR_RT, k)
            if TE_clade:
                LTR_RT_comp.append([k, v])
                if TE_clade.group("SUPFAM") == "Gypsy":
                    gypsy_count += 1
                elif TE_clade.group("SUPFAM") == "Copia":
                    copia_count += 1
                else:
                    unknown_count += 1
            else:
                non_LTR_RT_comp.append([k, v])

        # Write an overview regarding LTR-RTs in the EDTA annotation result
        with open(self._current_dir + "Module2_LTR_domain_classification/EDTA_LTR-RT_overview.txt",
                  mode='w') as overview:
            overview.write("\n\n%d LTR-RT sequence and %d non LTR-RT TE sequence parsed." % (
                len(LTR_RT_comp), len(non_LTR_RT_comp)))
            overview.write("Of all LTR-RT elements,\nGypsy: %d\nCopia: %d\nunknown: %d\n\n" % (
                gypsy_count, copia_count, unknown_count))

        # Output LTR-RT fasta and non LTR-RT fasta
        with open(self._current_dir + "Module2_LTR_domain_classification/EDTA_TElib_LTRRT.fa", mode='w') as output1:
            for indv_TE in LTR_RT_comp:
                output1.write(str(indv_TE[0]).rstrip("\n") + "\n")
                output1.write(str(indv_TE[1]).rstrip("\n") + "\n")
        with open(self._current_dir + "Module2_LTR_domain_classification/EDTA_TElib_nonLTRRT.fa", mode='w') as output2:
            for indv_TE in non_LTR_RT_comp:
                output2.write(str(indv_TE[0]).rstrip("\n") + "\n")
                output2.write(str(indv_TE[1]).rstrip("\n") + "\n")

    def TEsorter_deployer(self):
        print("Deploying TEsorter on separated LTR-RTs")
        os.system("cd %s && TEsorter %s -db rexdb-plant -p 20 && cd .." % (
            self._current_dir + "Module2_LTR_domain_classification",
            self._current_dir + "Module2_LTR_domain_classification/EDTA_TElib_LTRRT.fa"))

    def TEsorter_result_normalizer(self):
        print("Started parsing TEsorter cls.tsv file and generate csv files for further use.")
        TE_object_comp = []
        final_table_comp = []

        intact_LTR_crit = re.compile(r'^(?P<TE_ID>TE_\d+)_INT#LTR/.*?$')
        solo_LTR_crit = re.compile(r'^(?P<TE_ID>TE_\d+)_LTR#LTR/.*?$')
        other_TE_crit = re.compile(r'^(?P<TE_ID>TE_\d+)#.*?$')

        class TEObject(object):
            def __init__(self, tag, order, superfamily, clade, complete, strand, domains, domain_count,
                         LTR_identify_index,
                         TEID):
                self.tag = tag
                self.order = order
                self.superfamily = superfamily
                self.clade = clade
                self.complete = complete
                self.strand = strand
                self.domains = domains
                self.domain_count = domain_count
                self.LTR_identify_index = LTR_identify_index
                self.TEID = TEID

            def convert_completeness(self):
                """
                Only for LTR-RT identification.
                non LTR-RT elements will be shown as "Unknown"
                """
                if self.complete == "yes":
                    self.complete = 1
                elif self.complete == "no":
                    self.complete = 0

            def count_domains(self):
                if not self.domains == "none":
                    self.domain_count = len(str(self.domains).split("|")) - 1
                else:
                    self.domain_count = 0

            def get_TE_identity(self):
                if re.match(intact_LTR_crit, self.tag):
                    self.LTR_identify_index = 1
                    match_result = re.search(intact_LTR_crit, self.tag)
                    self.TEID = match_result.group("TE_ID")

                elif re.match(solo_LTR_crit, self.tag):
                    self.LTR_identify_index = 2
                    match_result = re.search(solo_LTR_crit, self.tag)
                    self.TEID = match_result.group("TE_ID")

                else:
                    self.LTR_identify_index = 0
                    match_result = re.search(other_TE_crit, self.tag)
                    self.TEID = match_result.group("TE_ID")

        def parse_TEsorter_tsv(TEsorter_tsv):  # Make result dataframes using original TEsorter cls.tsv file
            with open(TEsorter_tsv, mode='r') as lines:
                for line in lines:
                    if not line.startswith("#"):
                        line_data = line.split("\t")

                        if len(line_data) == 7:
                            TE_object_comp.append(
                                TEObject(line_data[0], line_data[1], line_data[2], line_data[3], line_data[4],
                                         line_data[5],
                                         line_data[6].rstrip("\n"), 0, 0, ""))
            for indv_TE in TE_object_comp:
                indv_TE.convert_completeness()
                indv_TE.count_domains()
                indv_TE.get_TE_identity()

                final_table_comp.append(
                    [indv_TE.tag, indv_TE.order, indv_TE.superfamily, indv_TE.clade, indv_TE.complete, indv_TE.strand,
                     indv_TE.domains, indv_TE.domain_count, indv_TE.LTR_identify_index, indv_TE.TEID])

            return pd.DataFrame(final_table_comp)

        final_df = parse_TEsorter_tsv(self._current_dir +
                                      "Module2_LTR_domain_classification/EDTA_TElib_LTRRT.fa.rexdb-plant.cls.tsv")
        final_df.to_csv(
            self._current_dir + "Module2_LTR_domain_classification/" + self._accession_name + "_TEsorter_final.csv")


mainAccession = IndividualAccession()
mainAccession.split_TElib()
mainAccession.TEsorter_deployer()
mainAccession.TEsorter_result_normalizer()
