---
name: shared-git-repo-race-hazard
description: ABC3b/ABC3cと共有しているのはworking treeだけでなく.git自体——他セッションのcommitで自分のHEADが動き、コミット直前のファイル内容が差し戻されることがある
metadata:
  type: feedback
---

2026-09-04、pGC Section3.lean のdocstring訂正コミット後、`git log -1 --stat`で
自分のコミットメッセージのはずが**別セッション(ABC3b/cのいずれか)の
コミットメッセージ**(「CorrHyp Track B: Spec K が...」)になっており、
diffには自分のSection3.lean/blocked-leaves.json変更が**そのセッションの
コミットに合流**していた。さらにその後、Section3.leanの内容が
`git status`上「作業ツリーで変更あり」と表示され、中身を見ると**自分が
直前に書いた訂正が消えて古い版に戻っていた**(HEADには正しい版が
コミット済みなのに、作業ツリーだけ古い版)。

**原因**: ABC3b/c とは working tree だけでなく **`.git`(index・HEAD参照）
そのものを共有している**——別セッションが `git commit` すると、その時点で
index に載っている全ファイル(自分が `git add` だけして未 commit だった
ものも含む)がそのセッションのコミットに巻き込まれる。さらに何らかの
操作(別セッションの `git checkout`/`reset` 等)で作業ツリーのファイルが
HEAD と無関係に上書きされることもある。

**How to apply**: (1) `git commit` 直後は必ず `git log -1 --stat` で
**自分の書いたメッセージ・自分のファイルだけ**が载っているか確認する。
違えば合流が起きている——内容自体が正しければ実害はないので、次の
コミットで訂正を続ければよい。(2) 何か訂正した直後にもう一度その
ファイルを触る前は `git status`/`git diff <file>` で作業ツリーが
HEAD と一致しているか確認する。ズレていたら(古い内容に差し戻って
いたら)`git restore <file>` で HEAD の内容に戻してから作業を続ける
——**中身を見ずに `git add`/`git commit` すると、直前の訂正を自分の
手で握り潰すことになる**。(3) `git add` は常に明示パスのみ
(既存の指示どおり)、かつ commit 後の検証を1ステップとして習慣化する。

関連: [[sibling-session-coordination-via-listagents]]
