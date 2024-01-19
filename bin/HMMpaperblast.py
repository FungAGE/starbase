import argparse
import requests
from bs4 import BeautifulSoup
import csv
import time

# TODO: add taxonomy filter?


def paperblast_search(file_path):
    # PaperBLAST search URL
    paperblast_url = "https://papers.genomics.lbl.gov/cgi-bin/hmmSearch.cgi"

    # Reading the HMM file
    with open(file_path, "rb") as file:
        files = {"hmmfile": (file_path, file)}

        # Submitting a file upload to PaperBLAST
        response = requests.post(paperblast_url, files=files, data={"Go": "Go"})

        # Checking if the request was successful (status code 200)
        if response.status_code == 200:
            return response.text
        else:
            print(f"Error: {response.status_code}")
            return None


def wait_for_page_load(html_content):
    # Implement a mechanism to wait for the page to load
    # You may need to adjust the waiting time based on the actual page load time
    time.sleep(5)  # Adjust the sleep time as needed

    return html_content


def parse_paperblast_results(html_content):
    # Parsing the HTML content using BeautifulSoup
    soup = BeautifulSoup(html_content, "html.parser")

    # Extracting relevant information based on style attributes
    results = []
    for p_element in soup.find_all(
        "p", style=lambda x: x and "margin-top: 1em; margin-bottom: 0em;" in x
    ):
        species = p_element.find_next("i").text.strip()
        gene_id = p_element.find_next("a").text.strip()
        evidence = p_element.find_next(
            "span", style="font-family: sans-serif; font-size: smaller;"
        ).text.strip()
        ncbi_prot_link = p_element.find("a")["href"]

        # Extracting additional information from the nested <ul> elements
        ul_element = p_element.find_next(
            "ul", style="margin-top: 0em; margin-bottom: 0em;"
        )
        li_element = ul_element.find_next("li", recursive=False)
        paper_title = li_element.find("a").text.strip()
        authors = (
            li_element.find_next("small", recursive=False)
            .find_next("a", recursive=False)
            .get("title")
        )
        nested_ul = li_element.find_next("ul", recursive=False)
        paper_excerpt = nested_ul.text.strip()
        results.append(
            {
                "species": species,
                "gene_id": gene_id,
                "evidence": evidence,
                "ncbi_prot_link": ncbi_prot_link,
                "authors": authors,
                "paper_title": paper_title,
                "paper_excerpt": paper_excerpt,
            }
        )

    return results


def save_results_to_tsv(results, output_file):
    with open(output_file, "w", newline="", encoding="utf-8") as tsvfile:
        fieldnames = [
            "species",
            "gene_id",
            "evidence",
            "ncbi_prot_link",
            "paper_title",
            "authors",
            "paper_excerpt",
        ]
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(results)


def main():
    parser = argparse.ArgumentParser(
        description="Submit a search to PaperBLAST and save results in a TSV file."
    )
    parser.add_argument("hmm_file", help="Path to the HMM file")
    parser.add_argument("output_file", help="Path to the output TSV file")

    args = parser.parse_args()

    # Submitting the HMM file to PaperBLAST
    html_content = paperblast_search(args.hmm_file)

    if html_content:
        # Waiting for the page to load
        html_content = wait_for_page_load(html_content)

        # Parsing the HTML content
        search_results = parse_paperblast_results(html_content)

        # Saving results to a TSV file with additional_info items as separate rows
        save_results_to_tsv(search_results, args.output_file)
        print(f"Results saved to {args.output_file}")
    else:
        print("Search failed.")


if __name__ == "__main__":
    main()
