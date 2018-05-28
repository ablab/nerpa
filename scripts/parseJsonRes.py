import json

def parseJson(filename):
    content = ""
    with open(filename, 'r') as f:
        content = f.read()
    parsejson = json.loads(content)
    if (len(parsejson["garlic_results"][0]["query_abstraction"]["node_strings"]) == 0):
        return
    print(len(parsejson["garlic_results"][0]["subject_abstraction"]["node_strings"][0]["chemical_nodes"]))
    score = parsejson["garlic_results"][0]["scores"]["relative_score"]
    match_cnt = len(parsejson["garlic_results"][0]["query_abstraction"]["node_strings"][0]["chemical_nodes"])
    elem_cnt = len(parsejson["garlic_results"][0]["subject_abstraction"]["node_strings"][0]["chemical_nodes"])
    return score, match_cnt, elem_cnt


print(parseJson("/home/olga/bio/NRP/data/testGR/BOW79/out/(2)_garlic.json"))
