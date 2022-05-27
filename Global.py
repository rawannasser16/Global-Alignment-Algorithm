f=open('blosum.txt','r')
x=f.readlines()
f.close()
l=[]  #TWO DEMINSION
h=[""] #READ ROW

for i in x:
   singleRow = i.split(" ") 
   for j in singleRow:
       if j !="" and j !="\n":
           h.append(j)
   l.append(h)
   h=[]
def Global_Aligment_Protein(first_seq_pro, second_seq_pro):
    len_first_seq = len(first_seq_pro)
    len_second_seq = len(second_seq_pro)
    max_value = 0
    maxcell = (0, 0)
    Aligment = [[0 for first_seq_pro in range(len_first_seq + 1)] for second_seq_pro in range(len_second_seq + 1)]
    for i in range(len_second_seq + 1):
        Aligment[i][0] = -i
    for j in range(len_first_seq + 1):
        Aligment[0][j] = -j
    for i in range(1, len_second_seq + 1, 1):
        for j in range(1, len_first_seq + 1, 1):
           
                First_aligment_index = l[0].index(first_seq_pro[j-1])
                Second_aligment_index=l[0].index(second_seq_pro[i-1])
                score =l[First_aligment_index][Second_aligment_index] 
                Aligment[i][j] = max(Aligment[i - 1][j] - 10, Aligment[i][j - 1] - 10, Aligment[i - 1][j - 1] + int(score))
                if max_value <= Aligment[i][j]:
                    max_value = Aligment[i][j]
                    maxcell = (i, j)

    print("This is Maximum value is", str(max_value))

    print("This is the Poistion of Maximum value in", str(maxcell))

    Gap_First_seq = ""
    Gap_Second_seq = ""
    Match = ""
    i=len_second_seq 
    j=len_first_seq
    while i>0 and j > 0:
        up = Aligment[i - 1][j] - 10
        left = Aligment[i][j - 1] - 10
        
        First_aligment_index = l[0].index(first_seq_pro[j-1])
        Second_aligment_index=l[0].index(second_seq_pro[i-1])
            
        value=l[First_aligment_index][ Second_aligment_index]
        print("(",str(first_seq_pro[j-1]),",",str(second_seq_pro[i-1]),")")
            
           
        print("the value ",str(value))
        score = l[First_aligment_index][ Second_aligment_index]
        
        diagonal = Aligment[i - 1][j - 1] +int(score)
        if Aligment[i][j] == diagonal:
            Gap_First_seq += first_seq_pro[j - 1]
            Gap_Second_seq += second_seq_pro[i - 1]
           
            if first_seq_pro[j-1]==second_seq_pro[i-1]:
                Match += "|"
            else:
                Match += " "
            i -= 1
            j -= 1
        elif Aligment[i][j] == up:
            Gap_First_seq += "-"
            Match += " "
            Gap_Second_seq += second_seq_pro[i - 1]
            i -= 1
        else:
            Gap_First_seq += first_seq_pro[j - 1]
            Gap_Second_seq += "-"
            Match += " "
            j -= 1

    
    while (i>0):
          Gap_First_seq+="-"
          Gap_Second_seq+=second_seq_pro[i-1]
          Match+=" "
          i-=1
    while(j>0):
        Gap_First_seq+=first_seq_pro[j-1]
        Gap_Second_seq+="-"
        Match+=" "
        j-=1

    Gap_First_seq = Gap_First_seq[::-1]
    Match = Match[::-1]
    Gap_Second_seq = Gap_Second_seq[::-1]

    print (Gap_First_seq)
    print(Match)
    print(Gap_Second_seq)
#-------------------DNA Global Aligment----------------------
def Global_Aligment(first_seq, second_seq):
    len_first_seq=len(first_seq)
    len_second_seq=len(second_seq)
    Aligment = [[0 for first_seq in range(len_first_seq+1)] for second_seq in range(len_second_seq+1)] 
    for i in range (len_second_seq+1):
        Aligment[i][0]=-i
    for j in range (len_first_seq+1):
        Aligment[0][j] =-j
    for i in range(1,len_second_seq+1,1):
        for j in range(1,len_first_seq+1,1):
            if first_seq[j-1] == second_seq[i-1]:
                score = 1
            else :
                score =-2
            Aligment[i][j] = max(Aligment[i-1][j]-1, Aligment[i][j-1]-1,Aligment[i-1][j-1]+score)
    print("The Value of Last Cell ",str(Aligment[len_second_seq][len_first_seq]))

    Gap_First =""
    Gap_Second=""
    Match=""
    i = len_second_seq
    j = len_first_seq
    while i > 0 and j > 0:
        up = Aligment[i-1][j]-1
        left=Aligment[i][j-1]-1
        if  first_seq[j-1]==second_seq[i-1]:
            score=1
        else :
            score=-2
        diagonal = Aligment[i-1][j-1]+score
        if Aligment[i][j]==diagonal:
            Gap_First+=first_seq[j-1]
            Gap_Second+=second_seq[i-1]
            if score==1:
                Match+="|"
            else:
                Match+=" "
            i-=1
            j-=1
        elif Aligment[i][j] == up:
            Gap_First+="-"
            Match+=" "
            Gap_Second+=second_seq[i-1]
            i-=1
        else:
            Gap_First+=first_seq[j-1]
            Gap_Second+="-"
            Match+=" "
            j-=1

    while (i>0):
        Gap_First+="-"
        Gap_Second+=second_seq[i-1]
        Match+=" "
        i-=1
    while(j>0):
        Gap_First+=first_seq[j-1]
        Gap_Second+="-"
        Match+=" "
        j-=1
        
    Gap_First=Gap_First[::-1]
    Match=Match[::-1]
    Gap_Second=Gap_Second[::-1]
    
    print (Gap_First)
    print(Match)
    print(Gap_Second)

# calling the function

choice=""
first_seq=""
second_seq =""
first_seq_pro=""
second_seq_pro=""
print("Your choice may be DNA or Protein")
choice=input("Enter choice")
if choice=="DNA" or choice=="dna":
   print("Enter the First Sequence of DNA") 
   first_seq=input()
   print("Enter the Second Sequence of DNA") 
   second_seq=input()
   Global_Aligment(first_seq, second_seq)
elif choice=="Protein" or choice=="protein":
     print("Enter the First Sequence of protein") 
     first_seq_pro=input()
     print("Enter the Second Sequence of protein") 
     second_seq_pro=input()
     Global_Aligment_Protein(first_seq_pro, second_seq_pro)
    
else:
 print("Invalid choice")

