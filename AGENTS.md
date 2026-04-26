## Git policy

I prefer to perform all Git operations manually.

You may use Git only for read-only inspection. Permitted examples include:

- `git status`
- `git diff`
- `git diff --staged`
- `git log`
- `git show`

Do not run any Git command that modifies repository state, including the working tree, staging area, commits, branches, tags, remotes, refs, Git config, stashes, submodules, or worktrees.

Do not create, amend, rewrite, delete, merge, rebase, cherry-pick, revert, tag, fetch, pull, push, stash, stage, unstage, switch branches, check out files, reset files, clean files, or otherwise mutate Git state.

Do not run commands including, but not limited to:

- `git add`
- `git restore`
- `git checkout`
- `git switch`
- `git reset`
- `git clean`
- `git stash`
- `git commit`
- `git commit --amend`
- `git rebase`
- `git merge`
- `git cherry-pick`
- `git revert`
- `git tag`
- `git branch`
- `git fetch`
- `git pull`
- `git push`
- `git worktree`
- `git submodule`
- `git config`

Do not change the staging area unless I explicitly ask you to do so. New and modified files should remain unstaged.

At the end of a task, summarize the files changed and suggest a commit message, but do not stage, commit, or push anything.
