# numc

### Provide answers to the following questions for each set of partners.
- How many hours did you spend on the following tasks?
  - Task 1 (Matrix functions in C): 2 hrs
  - Task 2 (Writing the Python-C interface): 8 hrs
  - Task 3 (Speeding up matrix operations): 24 hrs
- Was this project interesting? What was the most interesting aspect about it?
  - <b> Yes, this project is interesting not just in figuring the best method to speed up the program, but it is also interesting that we can write a Python-C interface to get a basic understanding of how to generate a python object. The most interesting part I found in this project is that switching one method or tweaking the codes a little bit can lead to a great performance in speed. This happened when we just put more unrollings on MatrixMul. The speed went from 70x to 110x. </b>
- What did you learn?
  - <b> I learned a lot in this project, not just in terms of knowing how to speed up the program but I learned a lot on how to discover bugs using a randomly generated matrix and a reference. I also learned that cache is the key to speeding up most of the functions. What surprised me the most was that parallelism is not always the best method to use. By having parallelism for a small matrix, it is running even slower. So I discovered using #pragma omp parallel for if(condition) would help disable parallelism if the matrix size is too small. I also learned that for MatrixPow, it is not necessary to have parallelism and better caching. This is because we found that using repeated square would result in the runtime to be log(N). The only thing that is useful for speed up is, perhaps, copying the temp_result from MatrixMul to the current result, so that it can pass on to the next iteration. I also learned that having recursion is also not a great way to approaches these problems. One reason is that it opens up lots of frames, some can lead to bad caching, and unable to perform parallelism. </b>
- Is there anything you would change?
  - <b>Yes, I realized that MatrixPow does not need two Mul_Matrix. What I was thinking was that I need one to record the current result and the other to record the cache of the square. If I did it backward, I might not need the cache of the square, I just need to skip one squared matrixed if the pow is divisible by 2. However, I think the performance of MatrixPow has reached and no need to further optimize it. </b>

### If you worked with a partner:
- In one short paragraph, describe your contribution(s) to the project.
  - <b>I am mainly in charge of writing task2 and implementing task3's Matrix's pow, but I did task1 with my own version. I was also in charge of researching effective algorithms and proof-checking. Besides these, I am also in charge of debugging major bugs and seeking minor bugs discovered from custom tests. I came up with the idea of using bit representation to do MatrixPow. This reduced the runtime to approximate log(N).</b>

- In one short paragraph, describe your partner's contribution(s) to the project.
  - <b>My partner, Luna, was in charge of implementing task1 and task3. She implemented all functions with data-level parallelism. When we are stuck for MatrixMul, she would go to office hours and ask for help. She came up with an idea to use cache blocking for MatrixMul. This is why our Matrix multiplication is so fast.</b>
