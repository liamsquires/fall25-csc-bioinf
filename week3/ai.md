# AI Usage Disclosure

## AI Tool Used
I used Claude Code CLI with Claude Sonnet 4.

## Project Context
The AI was provided with a CLAUDE.md file containing:
- Project overview and goals
- Detailed Codon programming language hints and patterns
- Specific instructions for Python-Codon bridge configuration
- File structure and command references

I've included the CLAUDE.md file in the repo for reference.

## User Prompts
The following prompts were provided to the AI during this assignment:

1. "I've made some edits to Claude.md to describe my goal. Please read it and let's make a plan to write the codon code one to pass one test at a time."

2. "Why are you removing these?" (regarding type annotations)

3. "Why are the files .py instead of .codon? Also can you move them to a folder called codon_source?"

4. "My mistake. The imports were correct and in accordance with the hints in CLAUDE.md. I have returned them to what they should be. Undo the changes you just made if you think they may cause issues down the line."

5. "Awesome. I believe the remaining issue is that codon is not set up in the venv, so it's using codon from outside the venv, where numpy doesn't exist. How might I resolve this?"

6. "Could this have to do with the hint about TreeNode hashing"

7. "Perfect. I now have a couple tasks to finish things off. 1. Can you add a timer to the test script and have it output the time it takes to run all tests? 2. Could you make a script, evaluate.sh, which runs the tests in python and codon, and outputs the time taken like so: Language Runtime --- python 2000ms codon 1000ms"

8. "Can you now create an ai.md that discloses use of AI for this assignment - state that claude code was used, it's version, mention CLAUDE.md, and list the prompts that I sent."
