# recursive permutations in python

from functools import reduce

def permute(in_list):
    #print('inner',in_list)
    if len(in_list)==0:
        return [[]];


    # need to wrap this in a reduce
    #input() 
    # using python filter:  will replace with right-fold in next version

    m1 = map(lambda x: map(lambda y: [x]+y,permute(list(filter(lambda z: z !=x,in_list)))), in_list)

    m1_list = list(map(lambda x: list(x),m1))
    #print('mi1_list',m1_list)

    #return m1_list
    
    # return tidy list
    return reduce(lambda a,b : a+b, m1_list)


# test
L=[1,2,3]

print(permute(L))
