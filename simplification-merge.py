#chrom /pos /type / taille/nombre lectures tot / nombre de lectures support / fréquence allélique

with open('simplified_combiSV.vcf', 'r') as input_file, open('Results_combiSV.simplified.vcf', 'w') as output_file:
    is_first_line = True
    for line in input_file:
        if not line.startswith("##"):
        # Diviser la ligne en colonnes en utilisant la tabulation comme séparateur
            columns = line.strip().split('\t')
            print(columns)
            
            if len(columns) >= 9:
            # Supprimer les colonnes qui ne nous intéresse pas
                del columns[2]
                del columns[2]
                del columns[2]
                del columns[2]
                del columns[2]
                del columns[3]


                # Extraire des éléments (type, taille, nombre lectures total, nombre de lectures support, fréquence allélique) de la colonne 3:  INFO
                if is_first_line:
                    # Si c'est la première ligne, conserver la colonne telle quelle
                    is_first_line = False
                else:
                    info = columns[2].split(";")
                    extracted_info = []
                    for item in info:
                        if item.startswith(("SVTYPE=")):
                            extracted_info.append(item)
                        if item.startswith(("SVLEN=")):
                            extracted_info.append(item)
                        if item.startswith(("SUPPORT=")):
                            extracted_info.append(item)
                        if item.startswith(("COVERAGE=")):
                            extracted_info.append(item)
                        if item.startswith(("AF=")):
                            extracted_info.append(item)
                        if item.startswith(("SVCALLERS=")):
                            extracted_info.append(item)

                
                        # Reconstruire la colonne numéro 3
                        columns[2] = ";".join(extracted_info)

                output_line = '\t'.join(columns) + '\n'
                output_file.write(output_line)


