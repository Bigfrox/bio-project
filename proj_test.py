
"""
PROJECT
2016253072 명수환
"""


from datetime import datetime
import sys
import math
import matplotlib.pyplot as plt
import time
from sklearn.metrics import roc_curve, roc_auc_score


def getDataFromFile(filename):
    file = open(filename, "r")

    ontology_MF = dict()
    ontology_BP = dict()

    line = None
    line = file.readline()
    # print(line)
    while line != "":
        # print(line)
        while line != "[Term]\n":
            line = file.readline()  # * Term 단위
            if line == "":
                return ontology_MF, ontology_BP

        line = file.readline().split()
        # print(line)
        # * Read line after "[Term]\n"
        if line[0] == "id:":
            id = line[1]
            # print("id:",id)
            line = file.readline().split()

        if line[0] == "name:":
            # print("name:",line[1:])
            line = file.readline().split()

        if line[0] == "namespace:":
            namespace = line[1]
            if namespace != "molecular_function" and namespace != "biological_process":
                line = file.readline().split()
                continue
            # print("namespace:",namespace)#* molecular_function,biological_process
            line = file.readline().split()

        is_a_relation = set()
        part_of_relation = set()
        is_obsleted = False
        while line:
            # print(line)
            if line[0] == "is_obsolete:":
                if line[1] == "true":
                    # print(id,"is obsolete")
                    line = file.readline().split()
                    is_obsleted = True
                    break
            elif line[0] == "is_a:":
                is_a_relation.add(line[1])
                line = file.readline().split()
            elif line[0] == "relationship:":
                if line[1] == "part_of":
                    part_of_relation.add(line[2])  # * id
                    line = file.readline().split()
                else:
                    line = file.readline().split()
            else:
                line = file.readline().split()

        if is_obsleted:
            continue

        error_relation = is_a_relation.intersection(part_of_relation)
        if error_relation:
            print(id, "has two relationship - ", error_relation)

        # * Classify - BP,MF
        if namespace == "molecular_function":
            ontology_MF[id] = list(
                is_a_relation.union(part_of_relation)
            )  # * {id:relation}
        elif namespace == "biological_process":
            ontology_BP[id] = list(is_a_relation.union(part_of_relation))

    return ontology_MF, ontology_BP


def getTermDepth(ontology, start_node, root_node):
    shortestPathLength = -1  # * 초기값
    if start_node == root_node:
        return 0
    for v in ontology[start_node]:
        edges_count = 0
        parent = v
        # print("parent : ", parent)
        edges_count += 1
        if parent != root_node:
            start_node = parent
            edges_count += getTermDepth(ontology, start_node, root_node)
        else:  # * parent is root
            return 1
        # print("short: ",shortestPathLength)
        # print("edge count : ",edges_count)
        if shortestPathLength != -1:
            if edges_count < shortestPathLength:
                shortestPathLength = edges_count
        elif shortestPathLength == -1:
            shortestPathLength = edges_count
        else:  # * shortestPathLength < edges_count
            pass

    # print("short: ",shortestPathLength)
    return shortestPathLength


def getDataFromGAF(filename, BP, MF):

    file = open(filename, "r")

    line = None

    while True:

        line = file.readline()

        if not line:
            break

        if line[0] == "!":
            continue

        line = line.split("\t")

        term = line[4].strip()
        gene = line[2].strip()
        code = line[6].strip()
        confirm = line[3].strip()
        ontology = line[8].strip()

        if gene == "":

            continue
        if term == "":

            continue
        if code == "":

            continue
        if code == "IEA":

            continue
        if confirm[:3] == "NOT":

            continue
        if ontology == "C":

            continue

        if ontology == "F":

            MF[term].add(gene)

        elif ontology == "P":
            BP[term].add(gene)

    return BP, MF
def GetPPI_with_SORA(filename, annotation_BP, annotation_MF, ontology_BP, ontology_MF, annotation_dict_BP,annotation_dict_MF):
    global root_BP
    global root_MF
    
    file = open(filename, "r")
    line = None
    gene_sim = list()
    line_count = 0
    total_line = 307679
    #* find the longest path in shortest path(depth)
    
    longest_path_length_BP = 0 # init value
    for v in ontology_BP:
        tmp = getTermDepth(ontology_BP,v,root_BP)
        if tmp >= longest_path_length_BP:
            
            longest_path_length_BP = tmp

    longest_path_length_MF = 0 # init value
    for v in ontology_MF:
        tmp = getTermDepth(ontology_MF,v,root_MF)
        if tmp >= longest_path_length_MF:
            
            longest_path_length_MF = tmp
    while True:
        line_count += 1
        
        if line_count % int(total_line / 1000) == 0:
            progress = round(line_count * 100 / total_line, 3)
            print("Progress : {0}%".format(progress))
            if progress >= 0.5:
                return gene_sim
            
        line = file.readline().split()
        print("line : ", line)
        print("line count : ", line_count)
        if not line:
            break
        gene1 = line[0]
        gene2 = line[1]

        
        isinBP = False #* init value
        isinMF = False
        
        # * gene1과 gene2가 BP에 속했는지 MF에 속했는지, 혹은 어디에도 속하지 않은 경우인지.
        if gene1 in annotation_dict_BP and gene2 in annotation_dict_BP:
            isinBP = True
        if gene1 in annotation_dict_MF and gene2 in annotation_dict_MF:
            isinMF = True
        if not isinBP and not isinMF:
            continue

        
        sim_BP_list = list()
        sim_MF_list = list()
        

        if isinBP:
            #* BP에만 속한 경우
            for v1bp in annotation_dict_BP[gene1]:
                C1_BP = v1bp
                for v2bp in annotation_dict_BP[gene2]:

                    C2_BP = v2bp
                    sim_BP = GetSimilarityUsingSORA(
                        C1_BP, C2_BP, ontology_BP, annotation_BP,root_BP,longest_path_length_BP
                    )
                    sim_BP_list.append(sim_BP)
                    
                #max_BP_list.append(max(sim_BP_list))
            avg_BP = sum(sim_BP_list)/len(sim_BP_list)
            

        if isinMF:
            #* MF에만 속한 경우
            
            #print("longest term : ", longest_term)
            #print("longest path length : ", longest_path_length)
            for v1mf in annotation_dict_MF[gene1]:
                C1_MF = v1mf
                for v2mf in annotation_dict_MF[gene2]:

                    C2_MF = v2mf
                    
                    
                    sim_MF = GetSimilarityUsingSORA(
                        C1_MF, C2_MF, ontology_MF, annotation_MF,root_MF,longest_path_length_MF
                    )

                    sim_MF_list.append(sim_MF)
                #max_MF_list.append(max(sim_MF_list))
            avg_MF = sum(sim_MF_list)/len(sim_MF_list)
            

            

        if sim_BP_list and sim_MF_list:
            # * BP와 MF 둘 다 속한 경우
            
            if avg_BP > avg_MF:
                sim = avg_BP
            else:
                sim = avg_MF

        elif sim_BP_list and not sim_MF_list:
            
            sim = avg_BP

        elif sim_MF_list and not sim_BP_list:
            sim = avg_MF
        else:
            #* Error case
            continue

        similarity = sim
        
        gene_sim.append([gene1, gene2, similarity])
        
    
    
    file.close()
    return gene_sim

def GetSimilarityUsingSORA(C1, C2, ontology, annotation,root_node,longest_path_length):
        
    # print("C1 : ", C1)
    # print("C2 : ", C2)
    IC_C1 = GetInformationContentSORA(C1,ontology,root_node=root_node,longest_term_depth=longest_path_length)
    IC_C2 = GetInformationContentSORA(C2,ontology,root_node=root_node,longest_term_depth=longest_path_length)
    IC_C0 = 0 #* init value
    common_ancestors = GetCommonAncestor(C1, C2, ontology)
    if not common_ancestors:
        return 0
    for ancestor in common_ancestors:
        if ancestor not in annotation:
            continue
        #print("ancestor : ", ancestor)
        tmp = GetInformationContentSORA(ancestor, ontology,root_node,longest_path_length)
        
        if IC_C0 <= tmp:
            IC_C0 = tmp
            
            C0 = ancestor
    # print("IC(C0) : ",IC_C0)
    # print("IC(C1) : ",IC_C1, "IC(C2) : ", IC_C2)
    similarity = 2 * IC_C0 / (IC_C1+IC_C2) # 이 값이 1.0 이하라는게 보장이 되지 않는다.
    if similarity > 1.0:
        similarity = 1.0
    #print("similarity : ", similarity)
    return similarity



def GetAllAncestor(nodeset, ontology):
    if not nodeset:
        return set()
    newset = set()

    for node in nodeset:

        if ontology[node]:
            for v in ontology[node]:
                newset.add(v)

    newset = newset.union(GetAllAncestor(newset, ontology))

    return newset


def GetCommonAncestor(term1, term2, ontology):
    
    set1 = {term1}
    set2 = {term2}
    ancestor1 = GetAllAncestor(set1, ontology)
    ancestor2 = GetAllAncestor(set2, ontology)
    
    if not ancestor1:
        return set()
    if not ancestor2:
        return set()
    common = ancestor1.intersection(ancestor2)

    
    return common



def inferred(ontology,annotation,total_length,root_node):
    print("[*]inferred start")

    while True:
        
        
        if len(annotation[root_node]) == total_length:
            satisfied = True
        else:
            satisfied = False
        

        if satisfied:
            break

        if not satisfied:
            
            for child in annotation:
                
                for parent in ontology[child]:
                    if parent in annotation:
                        
                        annotation[parent] = annotation[parent].union(annotation[child])
                        

            
    
    print("[*]inferred end")    


def GetTotalLengthOfGenes(annotation,root_node):
    global root_MF
    global root_BP

    gene_set_ = set()
    for term in annotation:
        gene_set_ = gene_set_.union(annotation[term])
        
        
    
    if root_node == root_BP:
        global BP_gene_set_length
        BP_gene_set_length = len(gene_set_)
        print("BP len : ",BP_gene_set_length)
    if root_node == root_MF:
        global MF_gene_set_length
        MF_gene_set_length = len(gene_set_)
        print("MF len : ", MF_gene_set_length)


def GetInformationContentSORA(term,ontology,root_node,longest_term_depth):#*SORA
    #* Specificity(term) * Coverage(term)
    #* specificity는 getTermdepth로
    #* Coverage는 후손 수로 계산
    specificity = getTermDepth(ontology,term,root_node) / longest_term_depth
    coverage = GetCoverage(term,ontology)
    # print("[term] : ",term)
    # print("specificity : ", specificity)
    # print("cov. : ", coverage)
    
    return specificity*coverage

def GetCoverage(term,ontology):
    # * 1 - ((log(term의 후손 수) + 1 ) / log(total terms 수))
    descendants = GetAllDescendants(term,ontology)
    # print("Descendants : ", descendants)
    # print("length : ", len(descendants))
    # print("ontology log.length : ", math.log(len(ontology)))
    return 1 - ((math.log(len(descendants)+1)) / math.log(len(ontology)))

def GetAllDescendants(term,ontology):
    # print("TERM : ", term)
    descendants_set = set()
    for child in ontology:
        if term in ontology[child]:
            descendants_set.add(child)
    # print("desc set : ",descendants_set)
    return descendants_set


def roc_curve(gene_sim,SE,OMSP):
    threshold = 0 #* line number
    TP_num = 0
    FP_num = 0
    FN_num = 0
    TN_num = 0 #* init value
    for element in gene_sim[:1]: #*initialize
        
        if len(element) == 4:
            element.append(0) #* init value
        
        if element[3] == "P":
            element[4] = "TP"
            #print(element)
            TP_num += 1
        else: #* element[3] == "N"
            element[4] = "FP"
            #print(element)
            FP_num += 1
    #print("==================Threshold")

    for element in gene_sim[1:]: #*initialize
        if len(element) == 4:
            element.append(0) #* init value

        if element[3] == "P":
            element[4] = "FN"
            #print(element)
            FN_num += 1
        else: #* element[3] == "N"
            element[4] = "TN"
            #print(element)
            TN_num += 1

    SE.append(TP_num / (TP_num+FN_num))
    OMSP.append(FP_num / (TN_num+FP_num))

    AUC_value = SE[threshold] * OMSP[threshold] #* first value
    
    threshold = 1
    while True:
        #print("\n\n\ ================\n\n")
        if threshold+1 > len(gene_sim):
            break
        
        element = gene_sim[threshold]
        if len(element) == 4:
            element.append(0) #* init value
        
        if element[3] == "P":
            #if element[4] == "FN":
            FN_num -= 1
            element[4] = "TP"
            TP_num += 1
        else: #* element[3] == "N"
            TN_num -= 1
            element[4] = "FP"
            FP_num += 1
    
        SE.append(TP_num / (TP_num+FN_num))
        OMSP.append(FP_num / (TN_num+FP_num))
        print("= = = = = ")
        print(SE[threshold])
        print(OMSP[threshold])
        print("SE[threshold] * OMSP[threshold] : ", SE[threshold] * OMSP[threshold])
        AUC_value += SE[threshold] * OMSP[threshold]
        threshold += 1
    
    return gene_sim, SE, OMSP, AUC_value

#def roc_curve(gene_sim,SE,OMSP):
    threshold = 0 #* line number
    # SE_ = list()
    # OMSP_ = list()
    
    while True:
        TP_num=0
        FP_num=0
        FN_num=0
        TN_num = 0 #* init value
        #print("\n\n\ ================\n\n")
        if threshold+1 > len(gene_sim):
            break
        for element in gene_sim[:threshold+1]:
            
            if len(element) == 4:
                element.append(0) #* init value
            
            if element[3] == "P":
                element[4] = "TP"
                #print(element)
                TP_num += 1
            else: #* element[3] == "N"
                element[4] = "FP"
                #print(element)
                FP_num += 1
        #print("==================Threshold")

        for element in gene_sim[threshold+1:]:
            if len(element) == 4:
                element.append(0) #* init value

            if element[3] == "P":
                element[4] = "FN"
                #print(element)
                FN_num += 1
            else: #* element[3] == "N"
                element[4] = "TN"
                #print(element)
                TN_num += 1
        
        SE.append(TP_num / (TP_num+FN_num))
        OMSP.append(FP_num / (TN_num+FP_num))
        threshold += 1
    
    return gene_sim, SE, OMSP



    
def plot_roc_curve(fper, tper, method):
    plt.plot(fper, tper, color='red', label='ROC')
    plt.plot([0, 1], [0, 1], color='green', linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic Curve(' + method +')')
    plt.legend()
    plt.show()

def main():

    if len(sys.argv) != 5:
        print("No input file.")
        print("<Usage> project.py go.obo goa_human.gaf biogrid_human_ppi_cln.txt human_protein_complexes.txt")
        return -1

    input_filename = sys.argv[1]
    input_filename2 = sys.argv[2]
    input_filename3 = sys.argv[3]
    input_filename4 = sys.argv[4]
    

    # input_filename = "go.obo"
    # input_filename2 = "goa_human.gaf"
    # input_filename3 = "biogrid_human_ppi_cln.txt"

    start_time = datetime.now()

    ontology_MF, ontology_BP = getDataFromFile(input_filename)

    print("length of MF :", len(ontology_MF))
    print("length of BP :", len(ontology_BP))

    global root_MF
    global root_BP

    for v1 in ontology_MF:

        if ontology_MF[v1] == []:
            root_MF = v1
            print("root node in MF : ", v1)

    for v2 in ontology_BP:

        if ontology_BP[v2] == []:
            root_BP = v2
            print("root node in BP : ", v2)

    error_count = 0
    for id_bp in ontology_BP:

        if id_bp == root_BP:
            continue
        for u1 in ontology_BP[id_bp][:]:

            if u1 in ontology_MF:
                # * ERROR -  Relation u and  v :  u is in MF, v is in BP
                error_count += 1
                ontology_BP[id_bp].remove(u1)
                # print(id_bp,"에서 ",u1,"를 삭제하였습니다.")

    for id_mf in ontology_MF:

        if id_mf == root_MF:
            continue
        for u2 in ontology_MF[id_mf][:]:
            if u2 in ontology_BP:
                # * ERROR - Relation u and  v  :  u is in BP, v is in MF
                error_count += 1
                ontology_MF[id_mf].remove(u2)

                # print(id_mf,"에서 ",u2,"를 삭제하였습니다.")
    print("error count : ", error_count)

    BP_annotation = ontology_BP.copy()
    MF_annotation = ontology_MF.copy()
    for v in BP_annotation:
        BP_annotation[v] = set()
    for v in MF_annotation:
        MF_annotation[v] = set()  

    BP_annotation, MF_annotation = getDataFromGAF(
        input_filename2, BP_annotation, MF_annotation
    )
    #BP_annotation={key:value for key, value in BP_annotation.items() if value != set()}
    #MF_annotation={key:value for key, value in MF_annotation.items() if value != set()}
    


    
    print("=========================================================================")
    print("length of BP_annotation : ", len(BP_annotation))
    print("length of MF_annotation : ", len(MF_annotation))
    global MF_gene_set_length
    global BP_gene_set_length
    
    GetTotalLengthOfGenes(MF_annotation,root_MF)
    GetTotalLengthOfGenes(BP_annotation,root_BP)
    
    print("MF_gene_set_length : ", MF_gene_set_length)
    print("BP_gene_set_length : ", BP_gene_set_length)

    print("BEFORE root_MF labels length ", len(MF_annotation[root_MF]))
    print("BEFORE root_BP labels length ", len(BP_annotation[root_BP]))


    BP_inferred_annotation = BP_annotation.copy()
    MF_inferred_annotation = MF_annotation.copy()
    
    # * inferred
    inferred(ontology_BP,BP_inferred_annotation,BP_gene_set_length,root_BP)
    inferred(ontology_MF,MF_inferred_annotation,MF_gene_set_length,root_MF)
    
    print("end inferred")

    # * show annotation

    print("AFTER root_MF labels length ", len(MF_inferred_annotation[root_MF]))
    print("AFTER root_BP labels length ", len(BP_inferred_annotation[root_BP]))
    BP_annotation={key:value for key, value in BP_annotation.items() if value != set()}
    MF_annotation={key:value for key, value in MF_annotation.items() if value != set()}

    annotation_dict_BP = dict()
    for gene in BP_inferred_annotation[root_BP]:
        for term in BP_annotation:
            if gene in BP_annotation[term]:
                if not gene in annotation_dict_BP: #처음인경우
                    annotation_dict_BP[gene] = set()
                annotation_dict_BP[gene].add(term)

    annotation_dict_MF = dict()
    for gene in MF_inferred_annotation[root_MF]:
        for term in MF_annotation:
            if gene in MF_annotation[term]:
                if not gene in annotation_dict_MF: #처음인경우
                    annotation_dict_MF[gene] = set()
                annotation_dict_MF[gene].add(term)
    #gene에 대해 연관된 term들 출력
    #for v in annotation_dict_MF:
    #    print(v, annotation_dict_MF[v],'\n')

    


    print("===========================================")
    print("PPI READ start")
    gene_sim = GetPPI_with_SORA(input_filename3, BP_annotation, MF_annotation, ontology_BP, ontology_MF,annotation_dict_BP,annotation_dict_MF)
    print("PPI READ end")
    complex_file = open(input_filename4,'r')
    complex_list = list()
    print("[*]complex data file open")
    while True:
        line = complex_file.readline().split()
        if not line:
            break
        complex_list.append(line)
    complex_file.close()
    print("[*]complex data file read end")

    for gene_pair in gene_sim:
        Found = False
        for complex_element in complex_list:
            if gene_pair[0] in complex_element and gene_pair[1] in complex_element:
                #print("Found!")
                #print("complex : ", complex_element)
                gene_pair.append("P")
                Found = True
                break
        if not Found:
            gene_pair.append("N")
    print("[*] Pos/Neg is assigned to PPI . . .")
    # * gene_sim을 유사도 내림차순으로 정렬
    gene_sim.sort(key=lambda x : x[2], reverse=True)
    SE = list()
    OMSP = list()
    gene_sim, SE, OMSP, AUC_value = roc_curve(gene_sim,SE,OMSP)

    plot_roc_curve(OMSP,SE,"New Method")
    print("AUC_value : ",AUC_value)
    #* 유사도 분포도 출력
    similarity = list()
    for v in gene_sim:
        similarity.append(v[2])
    
    print("[+] Time Elapsed : ", datetime.now() - start_time, "microseconds")
    plt.hist(similarity)
    plt.show()


if __name__ == "__main__":
    main()
