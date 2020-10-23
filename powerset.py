# powersets through strict recursion in Python


# recursive powerset function

def powerSet(input):
   

    # at tree leaf, return leaf
    if len(input)==0:
        return [[]];
    
    # if not at a leaf, trim and recurse
    
    # recursion is illustrated as follows:
    # S = {1,2,3}
    # S_trim = S without first element:
    # {(),(2),(3),(2,3)}
    # S_trim concatenated with first element:
    # {(1),(1,2),(1,3),(1,2,3)}
    # we keep the character sliced from front and concat it 
    # with result of recursion

    # use map to apply concatentation to all output from powerset

    leading = (input[0])
    new_input = input[1:len(input)]


    ps1 = list((powerSet(new_input)))
    # concatenate over powerset-ed set
    ps2 = map(lambda x: [leading]+x,ps1) 

    ps_list = list(map(lambda x: list(x),ps2))

    return ps1+ ps_list    

# test 
  
#input = [1,2,3,4];

#out=(powerSet(input))
#print(len(out))
#print(out)
