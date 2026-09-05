# メタ backlog —— 進め方の改善の台帳

`meta-optimizer` は spawn をまたいで記憶を持たない。ここが**その外部記憶**である。
数学トラックが依存グラフを外に置いたのと同じ構造(`ResearchPaper/orchestration.md` §2）。

```
起動 → この台帳を読む → 未着手の最上位を 1 件 → 測る → 直す → 結果を書き戻す → 止まる
```

**常駐プロセスは要らない。** 次に起動した agent が続きから始められる。

## 台帳の規約

- 1 件 = 1 節。**状態**は `未着手 / 測定済 / 実装済(提案) / 採用 / 却下`。
- **測定は必ず数字で書く。**「遅い気がする」は書かない。
- 1 回の起動で**直すのは 1 件だけ**(同時に複数変えると効果が測れない)。
- 実装は隔離 worktree。**自分でマージしない。** 採否は本体セッションが書き込む。
- 副作用の確認は毎回:`node tools/check.mjs --brief` の NG 件数(基準 13)と
  `node tools/graph.mjs` のノード数が変わらないこと。

---

## M1. 使い捨てスクリプト `_*.mjs` の集約

- **状態**: 測定中(2026-09-05 第 1 回に割当)
- **観測**: `tools/` 421 本のうち **399 本が `_*.mjs`**。2026-09-03 の実測で 397 本 / 10,528 行。
  未追跡のまま作業ツリーに溜まっているものもある。
- **仮説**: 台帳の 1 行を書くたびに使い捨てスクリプトを作っていた。`tools/ledger.mjs` が
  その一部を CLI に集約した(**過去のメタ改善の実例**)。残りも分類すれば減らせる。
- **効果**: —
- **採否**: —

## M2. ゲートの実コスト

- **状態**: 測定中(2026-09-05 第 1 回に割当)
- **観測**: `lake build ABC3` は約 6900 jobs で数分。出力の大半は `Replayed`(ログ再生)。
  MCP `lean_check` は 0.01〜1 秒。
- **仮説**: 時間の大半は再証明ではなく**トレース照合**。ならば変更ファイルの下流だけを
  `tools/graph.mjs --from <rel>` で列挙してビルドすれば十分条件になるかもしれない。
- **現在の運用**: 実装時は対象モジュール単体、`lake build ABC3` は節目で 1 回(2026-09-05 決定)。
- **効果**: —
- **採否**: —

## M3. ★`frontier.mjs` は import 依存しか見ていない(最重要)

- **状態**: 測定中(2026-09-05 第 1 回に割当)
- **観測(実例)**:
  - `Skeleton/PGC/Section1Cor13.lean` は `startable`(下流 27)と出る
  - 中身の `inertia_recoverable` は **Prop 1.2 に完全に帰着する**
    (`Found/PGC/InertiaTransport.lean::inertia_recoverable_of_prop12`、第 995)
  - Prop 1.2 は `Skeleton/PGC/Section1.lean` にあり、`Section1Cor13.lean` は**それを import していない**
  - → **スケジューラが「着手可能」と言うノードが、数学的には塞がっている**
- **仮説**: `.needs` の `.otherPaper` / `.derivation` に「import に現れない依存」が
  既に書かれている。拾えば `startable` の判定が正しくなる。
- **効果**: —
- **採否**: —

## M4. ★「不在」の測定が陳腐化する(今日 4 件の実例)

- **状態**: 測定中(2026-09-05 第 1 回に割当)
- **観測**: 2026-09-05 の 1 日で、「mathlib に X は無い」という記録が **4 件覆った**:

  | 記録 | 実際 |
  |---|---|
  | `Field (ULift α)` が無い(`Check/PGC/RefutationAttempts.lean`、2026-08-14) | `ULift.field` は**在る** |
  | 副有限の連続コホモロジーが無い(`Skeleton/PGC/Section1.lean` の `.needs`、2026-09-04) | `continuousCohomology` は**在る** |
  | `Ẑ` は mathlib に不在(`Found/PGC/AbsGalSurjections.lean`、2026-09-05) | `ProfiniteGrp.ProfiniteCompletion.completion` で**在る**(群として) |
  | `CompactSpace Gal(Ω/F)` が出せない(2026-09-05、`infer_instance` 失敗) | `FieldTheory/Galois/Profinite.lean:329` に**無名 instance で在る** |

- **なぜ高い**: 「不在」を誤ると**既にある数学を数百行書き直す**。今日の 4 件はいずれも
  「書き直しかけた」段階で見つかった。
- **仮説**: 主張と検索パターンを構造化して持ち、`.cache/mathlib-index.txt` に対して
  一括再測定できる。器は `ResearchPaper/lean-ecosystem.json` に既にある。
- **効果**: —
- **採否**: —

## M5. 索引が無名 instance と `@[to_additive]` 生成名を載せない(M4 の根本原因)

- **状態**: 未着手(2026-09-05 発見)
- **観測**: `.cache/mathlib-index.txt` は `tools/decl-index.mjs` のソース解析で作るため、
  `instance : CompactSpace Gal(K/k) := …`(名前なし)が**索引に入らない**。
  M4 の 4 件目はこれが原因で「無い」と誤判定された。名前で grep しても出ず、
  ソース直読で初めて見つかった。
- **仮説**: `decl-index.mjs` を無名 instance(型を鍵にする)と `to_additive` 生成名に
  対応させれば、`lean-search` の誤判定が構造的に減る。
- **効果**: —
- **採否**: —

## M6. セッションの実挙動が観測できない

- **状態**: 未着手(2026-09-05 発見)
- **観測**: `meta-optimizer` が見えるのは git 履歴・`tools/`・`lean-idioms.md` などの
  **成果物**だけ。「同じ検索を何回繰り返したか」「どのツール呼び出しが無駄だったか」という
  **最も濃い証拠は会話ログの中**にあり、渡していない。
- **仮説**: 無駄だった操作を機械可読な形で書き出す仕組み(あるいは転記の要約を渡す)を
  作れば、改善の質が上がる。
- **リスク**: 転記は巨大なので、そのまま渡すとコンテキストを食い潰す。要約の設計が要る。
- **効果**: —
- **採否**: —
