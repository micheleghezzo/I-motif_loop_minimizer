def sequence_decomposition(sequence, upper_limit, loop_nucleotide_type):
    first_run = ''
    first_loop = 0
    second_run = ''
    second_loop = 0
    third_run = ''
    third_loop = 0
    fourth_run = ''
    typing_error_list = []
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
        elif i == loop_nucleotide_type:
            if second_run == '':
                first_loop = first_loop + 1
            elif third_run == '':
                second_loop = second_loop + 1
            elif fourth_run == '':
                third_loop = third_loop + 1 
            else:
                typing_error_list.append('only four runs and three loops are accepted')
                break
        else:
            typing_error_list.append('only C (cytosines) and '+loop_nucleotide_type+' are accepted')
            break
    if typing_error_list == []:
        if first_run == '':
            typing_error_list.append('first run is missing')
        if first_loop == 0:
            typing_error_list.append('first loop is missing')
        elif first_loop > upper_limit:
            typing_error_list.append('first loop must not be longer than '+str(upper_limit)+' nucleotides')
        if second_run == '':
            typing_error_list.append('second run is missing')
        if second_loop == 0:
            typing_error_list.append('second loop is missing')
        elif second_loop > upper_limit:
            typing_error_list.append('second loop must not be longer than '+str(upper_limit)+' nucleotides')
        if third_run == '':
            typing_error_list.append('third run is missing')
        if third_loop == 0:
            typing_error_list.append('third loop is missing')
        elif third_loop > upper_limit:
            typing_error_list.append('third loop must not be longer than '+str(upper_limit)+' nucleotides')
        if fourth_run == '':
            typing_error_list.append('fourth run is missing')
        if first_run != third_run and first_run != '' and third_run != '':
            typing_error_list.append('first and third runs must be equeal')
        if second_run != fourth_run and second_run != '' and fourth_run != '':
            typing_error_list.append('second and fourth runs must be equal')
    return first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, typing_error_list   

def loop_as_string(loop_integer, loop_nucleotide_type):
    loop_string = ''
    for i in range(loop_integer):
        loop_string = loop_string + loop_nucleotide_type
    return loop_string

def check(first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, loop_nucleotide_type, P_list):
    sequence = first_run + loop_as_string(first_loop, loop_nucleotide_type) + second_run + loop_as_string(second_loop, loop_nucleotide_type) + third_run + loop_as_string(third_loop, loop_nucleotide_type) + fourth_run
    user_answer = input("Check the folding of 5'-"+sequence+"-3'\nDoes it fold into a monomeric I-motif using all the expected cytosines? yes/no\n").lower()
    while user_answer != 'yes' and user_answer != 'no':
        print('Please, type only yes or no!')
        user_answer = input('Does it fold into a monomeric I-motif using all the expected cytosines? yes/no\n').lower()
    P_list_reducted = eliminate(P_list, user_answer, first_loop, second_loop, third_loop)
    return user_answer, P_list_reducted, sequence

def eliminate(P_list, user_answer, first_loop, second_loop, third_loop):
    P_list_reducted = []
    if user_answer == 'yes':
        for i in P_list:
            if i[0] < first_loop or i[1] < second_loop or i[2] < third_loop:
                P_list_reducted.append(i)
    elif user_answer == 'no':
        for i in P_list:
            if i[0] > first_loop or i[1] > second_loop or i[2] > third_loop:
                P_list_reducted.append(i)
    return P_list_reducted    

def find_next(P_list, minimum):
    P_list_equal_first_third_loop = []
    P_list_different_first_third_loop = []
    p = 0
    for i in P_list:
        if i[0] == i[2]:
            P_list_equal_first_third_loop.append(i)
        else:
            P_list_different_first_third_loop.append(i)
    if P_list_equal_first_third_loop != []:
        for i in P_list_equal_first_third_loop:
            p_yes = len(P_list) - len(eliminate(P_list, 'yes', i[0], i[1], i[2]))
            p_no = len(P_list) - len(eliminate(P_list, 'no', i[0], i[1], i[2]))
            if p_yes <= p_no and p_yes > p:
                next_sequence = i
                p = p_yes
            elif p_no < p_yes and p_no > p:
                next_sequence = i
                p = p_no
    elif len(minimum) == 2 and [minimum[0][0]-1, minimum[0][1], minimum[0][2]] in P_list or len(minimum) == 2 and [minimum[0][0], minimum[0][1], minimum[0][2]-1] in P_list or len(minimum) == 2 and [minimum[1][0]-1, minimum[1][1], minimum[1][2]] in P_list or len(minimum) == 2 and [minimum[1][0], minimum[1][1], minimum[1][2]-1] in P_list:
        if [minimum[0][0]-1, minimum[0][1], minimum[0][2]] in P_list:
            next_sequence = [minimum[0][0]-1, minimum[0][1], minimum[0][2]]
        elif [minimum[0][0], minimum[0][1], minimum[0][2]-1] in P_list:
            next_sequence = [minimum[0][0], minimum[0][1], minimum[0][2]-1]
        elif [minimum[1][0]-1, minimum[1][1], minimum[1][2]] in P_list:
            next_sequence = [minimum[1][0]-1, minimum[1][1], minimum[1][2]]
        elif [minimum[1][0], minimum[1][1], minimum[1][2]-1] in P_list:
            next_sequence = [minimum[1][0], minimum[1][1], minimum[1][2]-1]
    elif len(minimum) == 1 and [minimum[0][0]-1, minimum[0][1], minimum[0][2]] in P_list or len(minimum) == 1 and [minimum[0][0], minimum[0][1], minimum[0][2]-1] in P_list:
        if [minimum[0][0]-1, minimum[0][1], minimum[0][2]] in P_list:
            next_sequence = [minimum[0][0]-1, minimum[0][1], minimum[0][2]]
        elif [minimum[0][0], minimum[0][1], minimum[0][2]-1] in P_list:
            next_sequence = [minimum[0][0], minimum[0][1], minimum[0][2]-1]
    elif P_list_different_first_third_loop != []:
        for i in P_list_different_first_third_loop:
            p_yes = len(P_list) - len(eliminate(P_list, 'yes', i[0], i[1], i[2]))
            p_no = len(P_list) - len(eliminate(P_list, 'no', i[0], i[1], i[2]))
            if p_yes >= p_no and p_yes > p:
                next_sequence = i
                p = p_yes
            elif p_no > p_yes and p_no > p:
                next_sequence = i
                p = p_no
    return next_sequence

def screening():     
    try:
        upper_limit = int(input("Input the upper limit to the loops lenght\n"))
    except:
        upper_limit = 0
    while upper_limit < 2:
        print('Please, type only integers > 1 ')
        try:
            upper_limit = int(input("Input the upper limit to the loops lenght\n"))
        except:
            upper_limit = 0
    loop_nucleotide_type = input("Input the nucleotide type of the loops\n").upper()
    while loop_nucleotide_type not in ['T', 'A', 'G', 'U', 'N']:
        print('Please, type only T (thymine), A (adenine), G (guanine), U (uracil), N (generic nucloetide)')
        loop_nucleotide_type = input("Input the nucleotide type of the loops\n").upper()
    P_list = []
    for first in range(1,upper_limit+1):
        for second in range(1,upper_limit+1):
            for third in range(1,upper_limit+1):
                P_list.append([first, second, third])
    sequence = input("Input the starting sequence from 5' to 3' end\n").upper()
    first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, typing_error_list = sequence_decomposition(sequence, upper_limit, loop_nucleotide_type)
    while typing_error_list != []:
        print('Please, check the sequence: ')
        for i in typing_error_list:
            print(i)
        sequence = input("Input the starting sequence from 5' to 3' end\n").upper()
        first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, typing_error_list = sequence_decomposition(sequence, upper_limit, loop_nucleotide_type)
    screening_summary = {}
    step = 0
    minimum = []
    user_answer = input('Does it fold into a monomeric I-motif using all the expected cytosines? yes/no\n').lower()
    while user_answer != 'yes' and user_answer != 'no':
        print('Please, type only yes or no!')
        user_answer = input('Does it fold into a monomeric I-motif using all the expected cytosines? yes/no\n').lower() 
    P_list_reducted = eliminate(P_list, user_answer, first_loop, second_loop, third_loop)
    screening_summary['step'+str(step)] = ["5'-"+sequence+"-3'", user_answer, 'sequences in P_list from '+str(len(P_list))+' to '+str(len(P_list_reducted))] 
    P_list = P_list_reducted
    if user_answer == 'yes':
        minimum.append([first_loop, second_loop, third_loop])             
    while P_list != []:
        step = step + 1
        next_sequence = find_next(P_list, minimum)
        first_loop = next_sequence[0]
        second_loop = next_sequence[1]
        third_loop = next_sequence[2]
        user_answer, P_list_reducted, sequence = check(first_run, first_loop, second_run, second_loop, third_run, third_loop, fourth_run, loop_nucleotide_type, P_list)
        screening_summary['step'+str(step)] = ["5'-"+sequence+"-3'", user_answer, 'sequences in P_list from '+str(len(P_list))+' to '+str(len(P_list_reducted))]   
        P_list = P_list_reducted
        if user_answer == 'yes':
            for i in minimum:
                if i[0] >= first_loop and i[1] >= second_loop and i[2] >= third_loop:
                    minimum.remove([i[0], i[1], i[2]])
            minimum.append([first_loop, second_loop, third_loop])
    if minimum == []:
        print('There is no folding for this system')
    elif len(minimum) == 1:
        print('The minimum pattern for this system is: ')
        print("5'-"+first_run + loop_as_string(minimum[0][0], loop_nucleotide_type) + second_run + loop_as_string(minimum[0][1], loop_nucleotide_type) + third_run + loop_as_string(minimum[0][2], loop_nucleotide_type) + fourth_run+"-3'")
    else:
        print('The minimum patterns for this system are: ')
        for i in minimum:
            print("5'-"+first_run + loop_as_string(i[0], loop_nucleotide_type) + second_run + loop_as_string(i[1], loop_nucleotide_type) + third_run + loop_as_string(i[2], loop_nucleotide_type) + fourth_run+"-3'")
    print('Summary of the screening:')
    for i in screening_summary:
        print(i+': '+str(screening_summary[i]))
    return screening_summary 

screening_summary = screening()  