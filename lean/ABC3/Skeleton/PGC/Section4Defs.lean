import ABC3.Skeleton.PGC.Setup
import ABC3.Interface.PGC.LocalFieldData

/-!
# [pGC] §4 — Lemma 4.1 / Theorem 4.2 が語る対象の**定義**

主張の本体は `ABC3/Skeleton/PGC/Section4.lean` にある。
本ファイルはその**定義だけ**を持つ(`Section1Defs.lean`・`Section2Defs.lean` と同じ作法)。

## ★なぜ定義と主張を分けたか(2026-09-08)

`node tools/import-audit.mjs`(2026-09-07 実測)によると、`Skeleton/PGC/Section4.lean` は
定義 1 個(`RamificationFiltration.filt`)と定理 2 個が同居した本で、そのせいで
**`Found/PGC/` の 202 本すべて**が「import すると循環する」状態だった
(`Section4` は `Section3` → `Section2` → `Section1` に加えて
`Found/PGC/OpenSubgroupSpan` と `Found/PGC/RamificationNaturality` を引く)。

★★詰まっていた本当の原因は `FilteredGroup` が `Found/PGC/FilteredGroup.lean` に
**しか**無かったことで、それは 2026-09-08 に `Skeleton/PGC/Setup.lean` へ移した。
本ファイルはその上で定義だけを受け取る。★定義そのものは 1 文字も変えていない。

★★**このファイルは `Found/` を 1 本も import しない**——`Setup` と
`Interface/PGC/LocalFieldData` だけで立つ。
`Section3Defs.lean` の側は `IsUniformizing` が `‖·‖` を使う都合で
`Found/PGC/LocalFieldNorm` の `scoped instance` を必要とするが、こちらは不要である。

## ★完全修飾名について

`RamificationFiltration` は `ABC3.Interface.PGC` の型なので、`filt` は
**`namespace ABC3.Skeleton.PGC` の外側**で完全修飾名を使って書く必要がある
(内側で `def RamificationFiltration.filt` と書くと
`ABC3.Skeleton.PGC.RamificationFiltration.filt` になってしまう)。
`Section4.lean` に置いてあったときと同じ書き方であり、**完全修飾名は
`ABC3.Interface.PGC.RamificationFiltration.filt` のまま 1 文字も変わらない**。
-/

/-- `K` の絶対 Galois 群を、`RF` の与える高次分岐群のフィルトレーションで飾った
filtered group(`Skeleton/PGC/Setup.lean` の `FilteredGroup`)。

★`RamificationFiltration`(`Interface/`)自体は `FilteredGroup` を知らない
(`Interface/PGC/LocalFieldData.lean` 冒頭の規約により、`Interface` は
他所を import しない)。この橋渡しは両方を import できる `Skeleton/` に置く。 -/
noncomputable def ABC3.Interface.PGC.RamificationFiltration.filt {p : ℕ} [Fact p.Prime]
    (RF : ABC3.Interface.PGC.RamificationFiltration p) (K : ABC3.Skeleton.PGC.PAdicLocalField p) :
    ABC3.Skeleton.PGC.FilteredGroup :=
  { G := K.absGal, Gv := RF.Gv K, isClosed := RF.isClosed K,
    isNormal := RF.isNormal K, antitone := RF.antitone K }

/-- 台帳の付随宣言(橋渡しの `def` であり、原典の項目そのものではない)。 -/
def ABC3.Interface.PGC.RamificationFiltration.filt.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 4, item := "Section 2 (RamificationFiltration.filt)",
    sectionId := "setup-2-herbrand" }
