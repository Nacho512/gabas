
Before anything below: this repository's engineering conventions (naming, module architecture, the core.lua input-validation contract, dependency-direction rules, numerical rigor, testing discipline) live in PROJECT_CONVENTIONS.md at the repo root. Read that first -- it applies regardless of which AI tool is working on this code, and it explicitly covers the Lua side only (INTERNAL/LUA-PART-PROJECT/), not the LaTeX side. Everything below is specific to how Codex itself should work here.

# Codex-Specific Mandatory Working Rules
## These instructions are intended for Codex

For every single session from now on, ALWAYS FOLLOW THE NEXT INSTRUCTIONS, UNLESS INSTRUCTED OTHERWISE:

1. Do not ask for confirmation unless the decision is irreversible or changes public behavior. Make reasonable assumptions, state them briefly, implement, and test. ASK ME ONLY IF YOU HAVE ANY QUESTION REGARDING A CONFLICT BETWEEN THE CURRENT SITUATION AND THESE RULES.

2. Any new explanatory source-code comment authored by Codex must begin with 'Codex:' after the language’s comment delimiter.

3. If you can't fix a glitch after reasonable tries and time invested, and/or if it's a limit-fringe case of the algorithm itself somewhere else, just leave an honest note letting me know of the existence of the bug/glitch/irreparable error, so I may try to polish it at a later time.

4. Add only necessary deterministic regression tests in the project’s existing test location. Do not create a new test structure unless no suitable structure exists RUN THE TESTS.

5. PRESERVE BACKWARD COMPATIBILITY.

6. USE EXISTING PROJECT PATTERNS AND DEPENDENCIES

7. DO NOT KNOWINGLY BREAK EXISTING BEHAVIOR. IF COMPATIBILITY MUST CHANGE, ASK FIRST.

8. VERIFICATION:

8.1. Run the smallest relevant tests/linter.
8.2. FIX ANY AND EVERY NEW FAILURE CAUSED BY YOUR CHANGES. PERIOD.

9. WORKING STYLE:

Return only: root cause, concise summary, files changed, test results, and remaining risks.

10. DO NOT NARRATE ROUTINE EXPLORATION

11. BE CONCISE; DO NOT NARRATE ROUTINE FILE INSPECTION.

13. Whenever I ask you something, or when reporting back after an inspection or a fixing session, ANSWER IN AT MOST 300 WORDS; NO EXCEPTIONS.

14. NEVER PASTE ENTIRE FILES IN THE RESPONSE; EDIT THE FILES DIRECTLY AND SUMMARIZE THE PATCH.

15. Run the narrowest relevant test first. If it passes, run the affected package's test suite. RUN THE FULL REPOSITORY SUITE ONLY IF REQUIRED BY PROJECT INSTRUCTIONS OR IF THE CHANGE HAS BROAD IMPACT.

16. Follow the existing project pattern. NEVER INTRODUCE A NEW ABSTRACTION UNLESS THE CURRENT DESIGN CANNOT SATISFY THE ACCEPTANCE CRITERIA.

17. Make the smallest correct patch.

18. DO NOT REFACTOR UNRELATED CODE; NEVER.

19. Do not update dependencies.

20. DO NOT REFORMAT UNTOUCHED FILES.

21. You are invited to leave a comment EXLUSUVELY when you consider it strictly necessary, FOLLOWING RULE #2. If you consider that because behavior or setup changes merit it, CREATE OTHER KIND OF DOCUMENTATION ONLY AFTER ASKING ME AND GETTING MY AUTHORIZATION; NO EXCEPTIONS.

22. DO NOT REPLACE EXISTING ARCHITECTURE MERELY BECAUSE ANOTHER DESIGN IS CLEANER.

23. Ask me only when required information is unavailable, instructions conflict, or a decision could cause irreversible or public-facing consequences.


B"H.
