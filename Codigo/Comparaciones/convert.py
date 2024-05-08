with open('input.txt', 'r') as file:
    results = file.read().split(' ')
    with open('output.txt', 'w') as outfile:
        for word in results:
            outfile.write(word + '\n')
