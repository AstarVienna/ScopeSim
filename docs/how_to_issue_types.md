# How to deal with issue types

- An ![Epic] _should_ have sub-issues (at some point) and serves as an overview.
- The sub-issues linked to an ![Epic] _may_ be either ![Task], ![Feature], ![Discussion] or ![Bug]. Any of these _may_ also exist as a standalone issue.
- An ![Epic] _must not_ have another ![Epic] as a sub-issue.
- ![Epic] level issues _should_ be kept "above" ![Feature] level issues.
- A ![Task] _must not_ have an ![Epic] or another ![Task] as a sub-issue, but _may_ itself be a sub-issue of either ![Epic] or ![Feature]. A ![Discussion] _may_ be a sub-issue of a ![Task] (if something needs to be clarified before the task in question can be completed, that is too large to discuss in the task itself, or is more abstract).
- If it turns out we really need another level in between ![Task] and ![Epic], we should introduce a `Batch` level (not to be confused with a "Patch").
- A ![Feature] _should not_ have another ![Feature] as a sub-issue. Usually, those will be either ![Task] or ![Discussion].

- A ![Bug] _should not_ have sub-issues.

## How to use labels now

- The `enhancement` label _should_ primarily be applied to PRs, however it _may_ also be applied to ![Epic] type issues. It is redundant on ![Feature] type issues and therefore _should not_ be applied there.
- The `bugfix` label _must_ only be applied to PRs. For issues that serve as bug reports, the issue type ![Bug] _must_ be used.
- A PR that is linked to a ![Bug] type issue (and will close it), _should_ have the `bugfix` label.

## How to deal with "sister PRs" in separate repos

If something needs a PR in both ScopeSim and the IRDB, best practice is to create two separate issues (usually of ![Task] or ![Feature] type) in either repo and then a parent (![Epic] or ![Feature]) issue in one of the repos (depending on where the change originates from, conceptually). Sub-issues can be linked across repos and this will make it easier to remember both of them. The individual PRs can then be linked to their respective sub-issues in their repo, but no one PR will be linked to the parent issues, since it cannot be closed by a single PR, because it spans two (or more) repos. This also reduces the amount of "This PR corresponds to X in Y." that is required to understand the context.

In any case, the corresponding PRs _must_ have identical branch names, if one is required to run the tests in the other. Usually, the IRDB PR will need the ScopeSim PR, so the latter should be created first, or at least have the identically named branch pushed already.

## Ideas

Perhaps do add Batch?

Maybe User Story? Or is that an Epic? chk usual difference

Label feedback?

[Task]: https://img.shields.io/badge/Task-darkgoldenrod
[Bug]: https://img.shields.io/badge/Bug-firebrick
[Feature]: https://img.shields.io/badge/Feature-mediumblue
[Epic]: https://img.shields.io/badge/Epic-forestgreen
[Discussion]: https://img.shields.io/badge/Discussion-indigo
