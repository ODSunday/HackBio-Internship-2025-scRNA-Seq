# Stage 0 Capstone Project
This goal of this project is to enable interns:
- Develop essential technical writing skills, focusing on clear and concise communication of complex scientific concepts. Interns will learn how to document research findings, write reports, and create educational content that is accessible to both scientific and general audiences.
- Write a simple Python script for printing the names, slack username, country, 1 hobby, affiliations of people on your team and the DNA sequence of the genes they love most.

The selected topic for the article (technical) writing: **What single-cell data is teaching us about cancer evolution**, with the focus of connecting scRNA-seq to tumour heterogeneity and therapy resistance.

### SECTION A: PYTHON SCRIPT

member_team_asparagine = {
    "Ayşenur Akcan": {
        "Slack Username": "@Ayşenur Akcan",
        "Country": "Türkiye",
        "Hobby": "Mending something",
        "Affiliation": "Istanbul Technical University",
        "DNA Sequence of Favourite Gene": "TP3- ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACTACTTCCTGAAAACAACGTTCT"
    },
    "Oni David Sunday": {
        "Slack Username": "@Sunday",
        "Country": "Nigeria",
        "Hobby": "Reading",
        "Affiliation": "Wageningen University and Research",
        "DNA Sequence of Favourite Gene": "AKT1 - TCTTCTCCCACTCCCCAGCAAAGTCCCCCTTTTGTGAGTGTAGCTGCCAGTACCTAGGTGAATGGTTGACTCCCCTCGGAGCCTTCCT..."
    },
    "Noor Ul Ain Amir": {
        "Slack Username": "@Noor Ul Ain Amir",
        "Country": "Pakistan",
        "Hobby": "Baking",
        "Affiliation": "Forman Christian College, University",
        "DNA Sequence of Favourite Gene": "SHH - ACAAGCTCTCCAGGCTTGCTACCATTTAAAATCAGACTCTTTTTGTCTTTTGATTGCTGTCTCGCGACCCAACTCCGATGTGTTCCGTTATCAGCGGCCGGCAGCCTGCCATTCCAGCCCCT..."
    },
    "Ozlem Kalkan": {
        "Slack Username": "@anyavala...",
        "Country": "Ukraine- Turkey",
        "Hobby": "doing sport",
        "Affiliation": "Bonn University",
        "DNA Sequence of Favourite Gene": "8OHY - ATGGCTCAAAGCACTGTTTTGCCTATGCACTGTTTATACGGGATATTCCTCGAAGGCAACCTAAAGATTCAAAAGAATGATCAGGAAGGACTAAAGAAGTTTAAAGATAATATCAAAAAATT…"
    },
    "Daria Kriuchkova": {
        "Slack Username": "@Daria",
        "Country": "Ukraine",
        "Hobby": "Knitting",
        "Affiliation": "Kharkiv Polytechnic Institute",
        "DNA Sequence of Favourite Gene": "BRCA1 - ATGGAAGTTGTCATTTTGTGTTTCCAGGATTTATTTGCTCTTCGTGTCTTTGGGTAGCTGG…"
    },
    "Minenhle Mayisela": {
        "Slack Username": "@Minenhle Mayisela",
        "Country": "Eswatini and South Africa",
        "Hobby": "Reading",
        "Affiliation": "University of the Witwatersrand",
        "DNA Sequence of Favourite Gene": "LMNA - GACAAATTCCTTGACCCGAGGAGGATAGGGATGTGGCCTTCGGTCTTTCCTCGCAGCTCCGGGGCAAGCTAGGAGTGGGATGGAAGTCGAGGTCCCTAATTTTTTAAGGGGAGGGTGCGGGGAGAAGGGGTAGTATGCGGAAACAGAGCGGGTATGAAGCTGGCTAACGCCGCGCGCCCCCTCCCAGGACCCGCTCCTGCCCCGCGCCGG"

    }
}      

def team_members(file):
    for member in file.keys():
        print(member)
        
def team_member_username(file):
    for member, info in file.items():
        print(f'Team member :{member} and his/her Slack username is: {info["Slack Username"]}')

def team_member_country(file):
    for member, info in file.items():
        print(f'Team member :{member} and his/her Country is: {info["Country"]}')

def team_member_hobby(file):
    for member, info in file.items():
        print(f'Team member :{member} and his/her hobby is: {info["Hobby"]}')
        
def team_member_affiliation(file):
    for member, info in file.items():
        print(f'Team member :{member} and his/her affiliation is: {info["Affiliation"]}')
        
def team_member_favseq(file):
    for member, info in file.items():
        print(f'Team member :{member} and his/her fav DNA seqeunce is: {info["DNA Sequence of Favourite Gene"]}')

if __name__ == "__main__":
    team_members(member_team_asparagine)
    print()  # just for spacing
    team_member_username(member_team_asparagine)
    print()
    team_member_country(member_team_asparagine)
    print()
    team_member_hobby(member_team_asparagine)
    print()
    team_member_affiliation(member_team_asparagine)
    print()
    team_member_favseq(member_team_asparagine)


### SECTION B: ARTICLE WRITING
