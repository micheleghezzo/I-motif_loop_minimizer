def sequence_decomposition(sequence):
    first_run = ''
    first_loop = 0
    second_run = ''
    second_loop = 0
    third_run = ''
    third_loop = 0
    fourth_run = ''
    error_list = []
    for i in sequence:   
        if i == 'C':
            if first_loop == 0:
                first_run = first_run + i
            elif second_loop == 0:
                second_run = second_run + i
            elif third_loop == 0:
                third_run = third_run + i
            else:
                fourth_run = fourth_run + i                
        elif i == 'T':
            if second_run == '':
                first_loop = first_loop + 1
            elif third_run == '':
                second_loop = second_loop + 1
            elif fourth_run == '':
                third_loop = third_loop + 1 
            else:
                error_list.append('only four runs and three loops are accepted')
                break
        else:
            error_list.append('only C (cytosines) and T (thymines) are accepted')
            break
    if error_list == []:
        if first_run == '':
            error_list.append('first run is missing')
        if first_loop == 0:
            error_list.append('first loop is missing')
        elif first_loop > 4:
            error_list.append('first loop must not be longer than 4 nucleotides')
        if second_run == '':
            error_list.append('second run is missing')
        if second_loop == 0:
            error_list.append('second loop is missing')
        elif second_loop > 4:
            error_list.append('second loop must not be longer than 4 nucleotides')
        if third_run == '':
            error_list.append('third run is missing')
        if third_loop == 0:
            error_list.append('third loop is missing')
        elif third_loop > 4:
            error_list.append('third loop must not be longer than 4 nucleotides')
        if fourth_run == '':
            error_list.append('fourth run is missing')
        if first_run != third_run and first_run != '' and third_run != '':
            error_list.append('first and third runs must be equeal')
        if second_run != fourth_run and second_run != '' and fourth_run != '':
            error_list.append('second and fourth runs must be equal')
    return first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, error_list   

def loop_as_string(loop_integer):
    loop_string = ''
    for i in range(loop_integer):
        loop_string = loop_string + 'T'
    return loop_string

def check(first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, step, screening_summary):
    step = step + 1
    sequence = first_run + loop_as_string(first_loop) + second_run + loop_as_string(second_loop) + third_run + loop_as_string(third_loop) + fourth_run
    folding = input("Check the folding of 5'-"+sequence+"-3'\nDoes it fold into a monomeric I-motif using all the expected cytosines? yes/no\n").lower()
    while folding != 'yes' and folding != 'no':
        print('Please, type only yes or no!')
        folding = input('Does it fold into a monomeric I-motif using all the expected cytosines? yes/no\n').lower()
    screening_summary['step'+str(step)] = ["5'-"+sequence+"-3'", folding]  
    return folding, step, screening_summary

def screening():
    sequence = input("Input the starting sequence from 5' to 3' end\n").upper()
    first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, error_list = sequence_decomposition(sequence)
    while error_list != []:
        print('Please, check the sequence: ')
        for i in error_list:
            print(i)
        sequence = input("Input the starting sequence from 5' to 3' end\n").upper()
        first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, error_list = sequence_decomposition(sequence)
    screening_summary = {}
    minimum = []    
    step = 0
    first_loop_limit = 1
    second_loop_limit = 1
    limit = [1, 1, 1]
    folding = input('Does it fold into a monomeric I-motif using all the expected cytosines? yes/no\n').lower()
    while folding != 'yes' and folding != 'no':
        print('Please, type only yes or no!')
        folding = input('Does it fold into a monomeric I-motif using all the expected cytosines? yes/no\n').lower()
    screening_summary['step'+str(step)] = ["5'-"+sequence+"-3'", folding]   
    while folding == 'no' and first_loop < 4 or folding == 'no' and second_loop < 4 or folding == 'no' and third_loop < 4:
        limit = [first_loop, second_loop, third_loop]
        if first_loop < 4 or third_loop < 4:
            if first_loop < 4:
                first_loop = first_loop + 1 
            if third_loop < 4:
                third_loop = third_loop + 1
            folding, step, screening_summary = check(first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, step, screening_summary)
            if folding == 'yes':                 
                first_loop_limit = first_loop
        elif second_loop < 4:
            second_loop = second_loop + 1  
            folding, step, screening_summary = check(first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, step, screening_summary)
            if folding == 'yes':                                   
                second_loop_limit = second_loop
    if folding == 'yes' and first_loop > 1 and second_loop >= second_loop_limit and third_loop > 1:
        while folding == 'yes' and second_loop > second_loop_limit:             
            second_loop = second_loop - 1 
            folding, step, screening_summary = check(first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, step, screening_summary)
            if folding == 'no':
                second_loop = second_loop + 1
        folding = 'yes'
        while folding == 'yes' and first_loop > 1:             
            first_loop = first_loop - 1
            folding, step, screening_summary = check(first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, step, screening_summary)
            if folding == 'yes':                 
                if first_loop == 1:
                    minimum.append([first_loop, second_loop, third_loop])
                    first_loop = first_loop_limit
                    break
            else:
                if first_loop < first_loop_limit:
                    minimum.append([first_loop + 1, second_loop, third_loop])
                    first_loop = first_loop_limit
                else:
                    first_loop = first_loop + 1
        folding = 'yes'
        while folding == 'yes' and third_loop > 1:
            third_loop = third_loop - 1
            folding, step, screening_summary = check(first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, step, screening_summary)
            if folding == 'yes':              
                if third_loop == 1:
                    minimum.append([first_loop, second_loop, third_loop])
            else:
                third_loop = third_loop + 1  
                minimum.append([first_loop, second_loop, third_loop])
        for i in minimum:
            if i[0] >= limit[0] and i[1] >= limit[1] and i[2] >= limit[2]:
                first_loop = i[0]
                second_loop = i[1]
                third_loop = i[2]
            else:
                first_loop = limit[0]
                second_loop = limit[1]
                third_loop = limit[2]
            if second_loop < 4 and first_loop != 1 or second_loop < 4 and third_loop != 1:
                second_loop = second_loop + 1
                second_loop_limit = second_loop
                if first_loop != 1:
                    first_loop = first_loop - 1
                if third_loop != 1:
                    third_loop = third_loop - 1 
                folding, step, screening_summary = check(first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, step, screening_summary)
                if folding == 'no':
                    while folding == 'no' and second_loop < 4:
                        second_loop = second_loop + 1         
                        folding, step, screening_summary = check(first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, step, screening_summary)
                        if folding == 'yes':                     
                            second_loop_limit = second_loop
                            first_loop_reduction_folding = 'yes'
                            while first_loop_reduction_folding == 'yes' and first_loop != 1: 
                                first_loop = first_loop - 1 
                                first_loop_reduction_folding, step, screening_summary = check(first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, step, screening_summary)
                                if first_loop_reduction_folding == 'no':
                                    first_loop = first_loop + 1
                            third_loop_reduction_folding = 'yes'
                            while third_loop_reduction_folding == 'yes' and third_loop != 1:                               
                                third_loop = third_loop - 1 
                                third_loop_reduction_folding, step, screening_summary = check(first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, step, screening_summary)
                                if third_loop_reduction_folding == 'no':
                                    third_loop = third_loop + 1
                            minimum.append([first_loop, second_loop, third_loop])
                    if folding == 'no' and first_loop < 4 or folding == 'no' and third_loop < 4:
                        while folding == 'no' and first_loop < 4:
                            first_loop = first_loop + 1 
                            folding, step, screening_summary = check(first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, step, screening_summary)
                            if folding == 'yes':                          
                                second_loop_reduction_folding = 'yes'
                                while second_loop_reduction_folding == 'yes' and second_loop != second_loop_limit: 
                                    second_loop = second_loop - 1 
                                    second_loop_reduction_folding, step, screening_summary = check(first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, step, screening_summary)
                                    if second_loop_reduction_folding == 'no':
                                        second_loop = second_loop + 1
                                minimum.append([first_loop, second_loop, third_loop])
                        folding = 'no' 
                        while folding == 'no' and third_loop < 4:
                            third_loop = third_loop + 1 
                            folding, step, screening_summary = check(first_run, first_loop_limit, second_run, second_loop, third_run, third_loop, fourth_run, step, screening_summary)
                            if folding == 'yes':  
                                second_loop_reduction_folding = 'yes'
                                while second_loop_reduction_folding == 'yes' and second_loop != second_loop_limit:                                  
                                    second_loop = second_loop - 1 
                                    second_loop_reduction_folding, step, screening_summary = check(first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, step, screening_summary)
                                    if second_loop_reduction_folding == 'no':
                                        second_loop = second_loop + 1
                                minimum.append([first_loop, second_loop, third_loop])
                else:                   
                    first_loop_reduction_folding = 'yes'
                    while first_loop_reduction_folding == 'yes' and first_loop != 1:               
                        first_loop = first_loop - 1 
                        first_loop_reduction_folding, step, screening_summary = check(first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, step, screening_summary)
                        if first_loop_reduction_folding == 'no':
                            first_loop = first_loop + 1
                    third_loop_reduction_folding = 'yes'
                    while third_loop_reduction_folding == 'yes' and third_loop != 1: 
                        third_loop = third_loop - 1 
                        third_loop_reduction_folding, step, screening_summary = check(first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, step, screening_summary)
                        if third_loop_reduction_folding == 'no':
                            third_loop = third_loop + 1
                    minimum.append([first_loop, second_loop, third_loop])
    else:
        minimum.append([first_loop, second_loop, third_loop])
    if len(minimum) == 0:
        print('There is not folding for this system')
    elif len(minimum) == 1:
        print('The minimum pattern for this sistem is: ')
        print("5'-"+first_run + loop_as_string(minimum[0][0]) + second_run + loop_as_string(minimum[0][1]) + third_run + loop_as_string(minimum[0][2]) + fourth_run+"-3'")
    else:
        print('The minimum patterns for this sistem are: ')
        for i in minimum:
            print("5'-"+first_run + loop_as_string(i[0]) + second_run + loop_as_string(i[1]) + third_run + loop_as_string(i[2]) + fourth_run+"-3'")
    print('Summary of the screening:')
    for i in screening_summary:
        print(i+': '+str(screening_summary[i])) 
    return screening_summary 

screening_summary = screening()
