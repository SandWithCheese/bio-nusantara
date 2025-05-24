import pandas as pd
from pygbif import occurrences as occ


def download_gbif_data():
    taxa = {
        "Mammalia": 359,
        "Aves": 212,
        "Amphibia": 131,
        "Crocodylia": 11493978,
        "Squamata": 11592253,
        "Testudines": 11418114,
        "Sphenodontia": 11569602,
    }

    reptile_limit = 1500
    default_limit = 6000
    per_page = 300

    for taxon_name, class_key in taxa.items():
        max_records = (
            reptile_limit
            if taxon_name in ["Crocodylia", "Squamata", "Testudines", "Sphenodontia"]
            else default_limit
        )
        print(
            f"Downloading {taxon_name} data (classKey: {class_key}, max {max_records})..."
        )

        all_results = []
        offset = 0

        while offset < max_records:
            try:
                data = occ.search(
                    country="ID",
                    hasCoordinate=True,
                    classKey=class_key,
                    limit=per_page,
                    offset=offset,
                )

                results = data.get("results", [])
                if not results:
                    break

                all_results.extend(results)
                print(f"Downloaded {len(results)} records at offset {offset}")

                offset += per_page

                if len(all_results) >= max_records:
                    break

            except Exception as e:
                print(f"Error downloading {taxon_name} at offset {offset}: {e}")
                break

        all_results = all_results[:max_records]

        if all_results:
            df = pd.DataFrame(all_results)
            filename = f"../data/raw/gbif_{taxon_name.lower()}_indonesia.csv"
            df.to_csv(filename, index=False)
            print(f"Saved {len(df)} total records to {filename}")
        else:
            print(f"No data found for {taxon_name}")

    print("Data download complete!")


def combine_reptilian_data():
    reptilian_orders = ["crocodylia", "squamata", "testudines", "sphenodontia"]
    combined_data = []

    for order in reptilian_orders:
        filename = f"../data/gbif_{order}_indonesia.csv"
        try:
            df = pd.read_csv(filename)
            df["reptilian_order"] = order.capitalize()
            combined_data.append(df)
            print(f"Added {len(df)} records from {order}")
        except FileNotFoundError:
            print(f"File not found: {filename}")
        except Exception as e:
            print(f"Error reading {filename}: {e}")

    if combined_data:
        reptiles_df = pd.concat(combined_data, ignore_index=True)
        reptiles_df.to_csv("../data/gbif_reptilia_combined_indonesia.csv", index=False)
        print(f"Combined reptilian dataset saved with {len(reptiles_df)} total records")
        return reptiles_df
    else:
        print("No reptilian data to combine")
        return None


if __name__ == "__main__":
    download_gbif_data()
    combine_reptilian_data()
