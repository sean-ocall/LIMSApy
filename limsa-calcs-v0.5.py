############################################################
#
#  Limsa Output Manipulation Tool for Reid Group
#
#  Author: Sean O'Callaghan spoc@unimelb.edu.au
#
#  This Revision (0.3): 18-Aug-2015
#
############################################################


import os, math

class Compound(object):
    def __init__(self, name, mass, conc):
        self.__name = name
        self.__conc = []
        self.__mass = float(mass)
        self.__conc.append(float(conc))
        self.__lipid_class = self.find_class(name)
        self.__ion = self.find_ion(name)

    def add_conc(self, conc):
        self.__conc.append(float(conc))

    def remove_conc(self, conc):
        if conc in self.__conc:
            self.__conc.remove(conc)

    def replace_concs(self, conc_list):
        self.__conc = conc_list

    def get_name(self):
        return self.__name

    def replace_name(self, name):
        self.__name = name

    def get_mass(self):
        return self.__mass

    def get_conc(self):
        return self.__conc

    def get_class(self):
        return self.__lipid_class

    def get_ion(self):
        return self.__ion

    def get_avg_conc(self):
        return sum(self.__conc)/len(self.__conc)

    def get_stdev_conc(self):
        sum_of_squares = 0.0
        avg_conc = self.get_avg_conc()
        for conc in self.__conc:
            result = conc - avg_conc
            result_squared = result*result
            sum_of_squares = sum_of_squares + result_squared

        stdev = math.sqrt(sum_of_squares/(len(self.__conc)-1))
        return stdev
     
    def find_class(self, name):
        if name[0:3] == 'Sph':
            lipid_class = 'SP'
            #print 'class is ', lipid_class
        elif name[0:3] == 'Cer':
            lipid_class = 'Cer'
        elif name[0:7] == 'Hex-Cer':
            lipid_class = 'Cer'
        elif name[0:6] == 'LacCer':
            lipid_class = 'Cer'
        elif name[0:10] == "(3'-sulfo)":
            lipid_class = 'Cer'
        elif name[0:2] == 'SM':
            lipid_class = 'SM'
        elif name[0:2] == 'MG':
            lipid_class = 'MG'
        elif name[0:2] == 'DG':
            lipid_class = 'DG'
        elif name[0:2] == 'TG':
            lipid_class = 'TG'
        elif name[0:3] == 'GPC':
            lipid_class = 'PC'
        elif name[0:3] == 'GPE':
            lipid_class = 'PE'
        elif name[0:3] == 'GPS':
            lipid_class = 'PS'
        elif name[0:3] == 'GPI':
            lipid_class = 'PI'
        elif name[0:3] == 'GPA':
            lipid_class = 'PA'
        elif name[0:3] == 'GPG':
            lipid_class = 'PG'
        elif name[0:4] == 'Chol':
            lipid_class = "Chol_and_CE"
        elif name[0:6] == 'OxChol':
            lipid_class = "Chol_and_CE"
        elif name[0:6] == "Cholic":
            lipid_class = "Bile_Acids"
        elif name[0:2] == "GT":
            lipid_class = "Gangliosides"
        elif name[0:2] == "GM":
            lipid_class = "Gangliosides"
        elif name[0:2] == "GD":
            lipid_class = "Gangliosides"
        elif name[0:2] == "GQ":
            lipid_class = "Gangliosides"
        elif name[0:3] == "Fuc":
            lipid_class = "Gangliosides"
        elif name[0:2] == "AC":
            lipid_class = "Acyl_Carnitines"
        elif name[0:4] == "OHAC":
            lipid_class = "Acyl_Carnitines"
        elif name[0:7] == "Retinol":
            lipid_class = "Retinoids"
        elif name[0:2] == "FA":
            lipid_class = "FA"
        elif name[0:4] == 'Acyl':
            lipid_class = "skin_SP"
        elif name[0:3] == 'GP1':
            lipid_class = 'GP'
        elif name[0:3] == 'GP2':
            lipid_class = 'GP'
        elif name[0:4] == 'GIPC':
            lipid_class = "plant_SP"
        elif name[0] == 'd':
            lipid_class = 'SP'
        elif name[0].lower() == 'x':
            lipid_class = 'product_ions'
        elif name[0:2] == 'CL':
            lipid_class = 'CL'            
        else:
            lipid_class = 'Unknown'
            print name, "not found"

        return lipid_class

    def find_ion(self, name):
        if 'NH4' in name:
            ion = 'NH4+'
        else:
            ion = 'MH+'

        return ion

                  
def limsa_calc():
    
    fileslist = []
    compounds = []

    for file in os.listdir("."):
        if file.endswith(".csv") and 'output' not in file and 'multiple' not in file:
            fileslist.append(file)
    print "filenames are:", fileslist

    for file in fileslist:
        fp = open(file, 'r')

        for i,line in enumerate(fp):
            parts = line.split(',')

            if i == 1:
                for j,part in enumerate(parts):
                    if part == 'Name':
                        name_idx = j
                    elif part == 'Mass':
                        mass_idx = j
                    elif part == 'Conc':
                        conc_idx = j
                    else:
                        pass

            if i>1:
                compound_idx = -1
                name = parts[name_idx]
                mass = parts[mass_idx]
                conc = parts[conc_idx]
                for k, compound in enumerate(compounds):
                    if compound.get_name() == name:
                        compound_idx = k

                if compound_idx == -1:
                    new_compound = Compound(name, mass, conc)
                    compounds.append(new_compound)
                else:
                    compounds[compound_idx].add_conc(conc)
        fp.close()

    return compounds



def correct_overlaps_2(compounds):
    
    mass_list = []
    for compound in compounds:
        mass = compound.get_mass()
        mass_list.append(mass)

    multiple_mass = []
    
    # Figure out the masses which appear multiple times
    for mass in mass_list:
        counter = mass_list.count(mass)
        if counter > 1 and mass not in multiple_mass:
            #print mass, "appears", counter, "times"
            multiple_mass.append(mass)
            
    #############################################################
    # Some code to print a useful intermediate file for debugging
    ############################################################
    op = open('multiples.csv','w')
    op.write('name, mass, conc1, conc2, conc3, conc4\n')
    for mass in multiple_mass:
        for compound in compounds:
            if compound.get_mass() == mass:
                op.write(compound.get_name() + ',' + str(mass) + ',')
                for conc in compound.get_conc():
                    op.write(str(conc) + ',')
                op.write('\n')
    op.close()
    ############################################################


    # If there are two or more compounds with the same mass
    for mass in multiple_mass:
        local_concs = []
        indices = []
        # Record the compound indices and concentrations
        for i, compound in enumerate(compounds):
            if compound.get_mass() == mass:
                local_concs.append(compound.get_conc())
                indices.append(i)
        
        # Let the first occurance of the mass take on
        # the concentration values and names of all occurances
        final_concs = local_concs[0]
        zero_concs = [0.0]
        for i, conc in enumerate(final_concs):
            if conc == 0:
                for concs in local_concs[1:]:
                    if concs[i] > 0:
                        final_concs[i] = concs[i]
                        zero_concs.append(0.0)

        new_name = ""
        for index in indices:
            # The new name is a concatenation of all of the names
            new_name = new_name + compounds[index].get_name() + "/"

        compounds[indices[0]].replace_concs(final_concs)
        compounds[indices[0]].replace_name(new_name)

        for index in indices[1:]:
            compounds[index].replace_concs(zero_concs)

def correct_for_overlaps(compounds):
    """
    LIMSA is producing output where a single mass value is randomly
    associated with two different lipid names. This causes gaps to
    appear in the concentration values and therefore zeros in the
    concentration columns in the output from this script.

    Here we try to fuse the two rows into a single row 

    """
    mass_list = []
    for compound in compounds:
        if compound.get_mass() not in mass_list:
            mass_list.append(compound.get_mass())

    for mass in mass_list:
        for i, compound in enumerate(compounds):
            if mass == compound.get_mass():
                concs = compound.get_conc()
                num_zeros = 0
                for conc in concs:
                    if conc == 0.0:
                        num_zeros = num_zeros + 1
                if num_zeros > 0:
                    for other_cpd in compounds[i:]:
                        if other_cpd.get_mass() == mass:
                            other_concs = other_cpd.get_conc()
                            non_zeros = []
                            for other_conc in other_concs:
                                if other_conc > 0:
                                    non_zeros.append(other_conc)

                            for non_zero in non_zeros:
                                compound.remove_conc(0.0)
                                compound.add_conc(non_zero)
                            

    

def lipid_class_calc(compounds, min_values=3):

    class_list = []
    for compound in compounds:
        if compound.get_class() not in class_list:
            class_list.append(compound.get_class())

    class_dict = {}
    for lipid_class in class_list:
        class_dict[lipid_class] = 0.0

    for lipid_class in class_list:
        for compound in compounds:
            test_for_zeros = sum(i>0.0 for i in compound.get_conc())
            if compound.get_class() == lipid_class and \
                    test_for_zeros >= min_values:
                class_dict[lipid_class] = class_dict[lipid_class] + \
                    compound.get_avg_conc()


    return class_dict

def write_output(compounds, outfilename, class_dict, min_values=3):
    
    op = open(outfilename, 'w')

    fileslist = []

    for file in os.listdir("."):
        if file.endswith(".csv") and 'output' not in file and 'multiple' not in file:
            fileslist.append(file.strip('.csv'))

    op.write('Output\n')
    op.write(',,,Conc, Conc, %Total, %Total, ')
    for name in fileslist:
        op.write('Conc,')
    op.write('\n')
    op.write('Name, Class, m/z, Mean, STDEV, Conc, STDEV,')
    for file in fileslist:
        op.write(file + ',')
    op.write('\n')

    class_list = []
    mass_list = []
    for compound in compounds:
        mass_list.append(compound.get_mass())
        if compound.get_class() not in class_list:
            class_list.append(compound.get_class())
    
    mass_list.sort()
    
    for lipid_class in class_list:
        # There is an issue with the same mass being recorded
        # for two completely different compounds
        # Cannot assume unique masses!
        # This would result in peaks being printed twice
        # if not caught by if statement
        for j,mass in enumerate(mass_list):
            for compound in compounds:
                if mass == compound.get_mass() and \
                        lipid_class == compound.get_class() and\
                        mass != mass_list[j-1]:
                
                    test_for_zeros = sum(i>0.0 for i in compound.get_conc())
                    if test_for_zeros >= min_values:
                        name = compound.get_name()
                        lipid_class = compound.get_class()
                        mass = compound.get_mass()
                        avg_conc = compound.get_avg_conc()
                        stdev_conc = compound.get_stdev_conc()
                        conc_list = compound.get_conc()
                        ion = compound.get_ion()

                        percent_total = 100*avg_conc/class_dict[lipid_class]
                        percent_stdev = 100*stdev_conc/class_dict[lipid_class]

                        op.write(str(name) + ',' + lipid_class + ','\
                                    + "%.4f"%mass + ',' + str(avg_conc) \
                                     + ',' + str(stdev_conc) \
                                     + ',' + str(percent_total) + ','+\
                                     str(percent_stdev) + ',')
                        for conc in conc_list:
                            op.write(str(conc) + ',')

                        if j+1 < len(mass_list) and \
                                mass_list[j] == mass_list[j+1]:
                            op.write('*,')
                            
                        op.write('\n')
        op.write(',,,,,Total for %s:,%s\n'%(lipid_class,class_dict[lipid_class]))

if __name__ == '__main__':
    compounds = limsa_calc()
    correct_overlaps_2(compounds)
    #correct_for_overlaps(compounds)
    dict_of_class_totals = lipid_class_calc(compounds, min_values=2)
    write_output(compounds, 'output.csv', dict_of_class_totals, min_values=2)
    
