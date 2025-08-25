
def WriteToOutputFile(pre_computed_text, file_name="output.txt"):
    try:
        with open(file_name, "w") as file:
            file.write(pre_computed_text)
        print(f"Text successfully written to {file_name}")
    except Exception as e:
        print(f"An error occurred: {e}")


