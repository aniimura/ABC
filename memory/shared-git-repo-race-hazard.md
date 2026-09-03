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

**続報(2026-09-04、`.lake/build`キャッシュでも類似の再現性トラブルを観測)**:
`hasSubst_g_subst_X`の一般化(型注釈なしでσを遅延推論に頼る既存3箇所)は
一度 `lake build ABC3`(6590 jobs)が完全成功した**直後**に、同じソースを
一切変更していないのに `lake build` を打ち直すと「typeclass instance
problem is stuck」で**再現性をもって**(2回連続で同一エラー)落ちる
現象に遭遇した。`git status`/`git diff`は該当ファイルにズレ無し
(HEADと一致)だったので `.git` 側の合流ではない——ABC3b/cと`lean/.lake/
build`(ビルドキャッシュ)も同じディレクトリを共有しているため、並行
セッションの `lake build` がoleanを書き換えている最中に自分のビルドが
古い/一部だけ更新された状態を読んだ可能性が疑われる(未確定)。
**How to apply**: 原因の特定に時間をかけるより、型注釈なしでの
`have h := generalizedLemma arg1 arg2`(戻り値の型から後方に多相型
引数を推論させる書き方)を型注釈つき(`have h : ExpectedType := ...`
または `generalizedLemma (σ := ConcreteType) arg1 arg2`)に直すほうが
安く恒久的——多相化した既存補題を呼ぶ箇所は型注釈を省略しない習慣に
切り替える。「一度通った」は「今後も通る」を保証しない(この環境では)。

関連: [[sibling-session-coordination-via-listagents]]・[[pgc-lubin-tate-existence-progress]]
