import os
import re
import shutil
import argparse

def clean_target_folder(target_folder):
    for filename in os.listdir(target_folder):
        if filename.endswith('.sav') or filename.endswith('.des.dat') or filename.endswith('.des.log'):
            file_path = os.path.join(target_folder, filename)
            os.remove(file_path)

def create_target_folder(source_folder, new_model_name):
    target_folder = os.path.join(os.path.dirname(source_folder), f"test_diode_{new_model_name}")
    os.makedirs(target_folder, exist_ok=True)
    return target_folder

def get_concentration_values(content):
    acceptor_match = re.search(r"\bAcceptor.*?Conc=([\d.e+-]+)", content, re.S)
    donor_match = re.search(r"\bDonor.*?Conc=([\d.e+-]+)", content, re.S)
    acceptor_conc = acceptor_match.group(1) if acceptor_match else None
    donor_conc = donor_match.group(1) if donor_match else None
    return acceptor_conc, donor_conc

def prompt_new_concentrations(old_acceptor, old_donor):
    new_acceptor = input(f"Enter new acceptor concentration (old: {old_acceptor}): ")
    new_donor = input(f"Enter new donor concentration (old: {old_donor}): ")
    return new_acceptor or old_acceptor, new_donor or old_donor

def replace_concentrations(content, old_acceptor, new_acceptor, old_donor, new_donor):
    # Replace the acceptor concentration specifically
    content = re.sub(rf"(Acceptor.*?Conc=){old_acceptor}", rf"Acceptor\g<1>{new_acceptor}", content, flags=re.S)
    # Replace the donor concentration specifically
    content = re.sub(rf"(Donor.*?Conc=){old_donor}", rf"Donor\g<1>{new_donor}", content, flags=re.S)
    return content

def get_last_voltage(content):
    voltage_matches = re.findall(r"Voltage\s*=\s*([-.\d]+)", content)
    return float(voltage_matches[-1]) if voltage_matches else None

def prompt_new_voltages(last_voltage):
    print(f"Last voltage in file: {last_voltage}")
    new_voltages = input("Enter additional voltages separated by commas (or press Enter to skip): ")
    return [float(v.strip()) for v in new_voltages.split(',')] if new_voltages else []

def add_voltage_blocks(content, new_voltages, model_name):
    # Remove the last '}' character to prevent out-of-scope errors
    content = content.rstrip()  # Strip any trailing whitespace
    if content.endswith('}'):
        content = content[:-1]  # Remove the last character if it's a closing brace

    for voltage in new_voltages:
        voltage_int = int(abs(voltage))  # Create suffix based on integer part
        new_block = f"""
    Save (FilePrefix="test_diode_{model_name}_{voltage_int}")
    Quasistationary 
    (Goal {{Name = "cathode" Voltage = {voltage}}} InitialStep = 1.0) {{
        Coupled (Iterations = 50 Digits = 5) {{Poisson electron hole}}
    }}
    Plot (FilePrefix="test_diode_{model_name}_{voltage_int}")
"""
        content += new_block

    content += '}\n'  # Add the closing brace at the end
    return content

def rename_and_copy_files(source_folder, target_folder, old_model_name, new_model_name):
    for filename in os.listdir(source_folder):
        if old_model_name in filename:
            new_filename = filename.replace(old_model_name, new_model_name)
            source_path = os.path.join(source_folder, filename)
            target_path = os.path.join(target_folder, new_filename)

            if filename == f"test_diode_{old_model_name}_des.cmd":
                with open(source_path, 'r') as file:
                    content = file.read()
                
                content = content.replace(old_model_name, new_model_name)
                old_acceptor, old_donor = get_concentration_values(content)
                if old_acceptor and old_donor:
                    new_acceptor, new_donor = prompt_new_concentrations(old_acceptor, old_donor)
                    content = replace_concentrations(content, old_acceptor, new_acceptor, old_donor, new_donor)
                
                last_voltage = get_last_voltage(content)
                if last_voltage is not None:
                    new_voltages = prompt_new_voltages(last_voltage)
                    if new_voltages:
                        content = add_voltage_blocks(content, new_voltages, new_model_name)

                with open(target_path, 'w') as file:
                    file.write(content)
            else:
                shutil.copy(source_path, target_path)

def main():
    parser = argparse.ArgumentParser(description="Process TCAD model files.")
    parser.add_argument("source", help="Path to the source folder with template files.")
    parser.add_argument("old_model_name", help="Old model name to replace.")
    parser.add_argument("new_model_name", help="New model name to use.")

    args = parser.parse_args()

    target_folder = create_target_folder(args.source, args.new_model_name)
    rename_and_copy_files(args.source, target_folder, args.old_model_name, args.new_model_name)
    clean_target_folder(target_folder)
    print(f"Files processed and saved to: {target_folder}")



if __name__ == "__main__":
    main()
