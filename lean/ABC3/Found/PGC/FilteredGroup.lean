import ABC3.Skeleton.PGC.Setup

/-!
# [pGC] Definition 2.3 — filtered group(**旧名の互換層**)

★★**2026-09-08: 中身は `ABC3/Skeleton/PGC/Setup.lean` へ移した。**
本ファイルに残っているのは「旧い完全修飾名 `ABC3.Found.PGC.FilteredGroup*` を
引き続き解決できるようにする `export`」だけである。

## なぜ移したか(測ってある)

`node tools/import-audit.mjs`(2026-09-07 実測)によると、

| Skeleton の本 | 定義 | 定理 | `Found/PGC` 202 本のうち import すると循環する本 |
|---|---|---|---|
| `PGC/Section1Defs` / `Section2Defs` / `Setup` | 5 / 2 / 5 | 0 | 0(0%) |
| `PGC/Section3`(`filteredGroupOf` `IsUniformizing`) | 2 | 2 | 200(99%) |
| `PGC/Section4`(`RamificationFiltration.filt`) | 1 | 2 | 202(100%) |

★詰まる原因は 1 つ、同じもの——**`FilteredGroup` が `Found` の本にしかない**こと。
`Section1Defs` / `Section2Defs` と同じ作法で `Section3Defs` / `Section4Defs` へ
定義を割っても、割った先が `ABC3.Found.PGC.FilteredGroup` を持ち越すので
**Skeleton の Defs が Found を import する**羽目になり、割る作業が 1 段で閉じない
(`import-audit --plan-all` が「★★これも毒 ⇒ この本も割る必要がある」と名指ししていた)。

移動先を `Interface/PGC/` ではなく `Skeleton/PGC/Setup.lean` にした理由
(G2 ゲートが witness を 2 つ要求すること・`Interface` の役割定義に合わないこと・
新しい import 辺が 1 本も増えないこと)は移動先の docstring に書いた。

## 本ファイルを消さなかった理由

* `import ABC3.Found.PGC.FilteredGroup` と書いている本が 3 つある
  (`Found/PGC/AbsGalRamificationFiltration.lean`・`Found/PGC/RamificationNaturality.lean`・
  `Skeleton/PGC/Section2.lean`)。消すと `lean/ABC3/Found.lean` も含めて
  書き換えが要る——`Found.lean` は他の agent が同時に編集している。
* `Found/PGC/AbsGalRamificationFiltration.lean::toFilteredGroup`(305 行目)は
  `namespace ABC3.Found.PGC` の内側で **`open ABC3.Skeleton.PGC` より前**に
  bare な `FilteredGroup` を書いている(その `open` は 409 行目)。
  下の `export` があればその本は 1 文字も直さずに通る。

★`export` は同じ定数への別名なので、`open ABC3.Skeleton.PGC ABC3.Found.PGC` を
同時に行っても曖昧にならない(`lean_check` 0.01 秒で確認済み)。
-/

namespace ABC3.Found.PGC

export ABC3.Skeleton.PGC (FilteredGroup FilteredGroup.Iso FilteredGroup.OuterIso)

end ABC3.Found.PGC
