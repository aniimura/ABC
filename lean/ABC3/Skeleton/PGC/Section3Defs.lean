import ABC3.Skeleton.PGC.Setup
import ABC3.Interface.PGC.LocalFieldData
import ABC3.Found.PGC.LocalFieldNorm

/-!
# [pGC] §3 — Corollary 3.1 / Definition 3.2 / Corollary 3.3 が語る対象の**定義**

主張の本体は `ABC3/Skeleton/PGC/Section3.lean` にある。
本ファイルはその**定義だけ**を持つ(`Section1Defs.lean`・`Section2Defs.lean` と同じ作法)。

## ★なぜ定義と主張を分けたか(2026-09-08)

`node tools/import-audit.mjs`(2026-09-07 実測)によると、`Skeleton/PGC/Section3.lean` は
定義 2 個と定理 2 個が同居した本で、そのせいで **`Found/PGC/` の 202 本のうち 200 本**が
「import すると循環する」状態だった。`Found` 側から `filteredGroupOf` /
`IsUniformizing` を引きたければ `Section3` ごと(= `Section2` → `Section1` →
`Found/PGC/ReciprocityDatumIndependence` …)引き込むしかなかったからである。

★★詰まっていた本当の原因は `FilteredGroup` が `Found/PGC/FilteredGroup.lean` に
**しか**無かったことで、それは 2026-09-08 に `Skeleton/PGC/Setup.lean` へ移した。
本ファイルはその上で定義だけを受け取る。★定義そのものは 1 文字も変えていない。

★名前空間は `Section3.lean` と同じ `ABC3.Skeleton.PGC` なので**完全修飾名は
1 文字も変わらない**(`Section1Defs` / `Section2Defs` のときと同じ)。

## ★`Found/PGC/LocalFieldNorm` を import している理由(実測)

`IsUniformizing` の本文は `‖x‖ = 1`(= `U_K`)を使う。`K.carrier` 上のノルムと
位相は `Found/PGC/LocalFieldNorm.lean` の **`scoped instance`**
(`normedField` / `normedAlgebra` …)が与えるものなので、その本を import して
`open ABC3.Found.PGC` しないと `IsUniformizing` は型が付かない:

```
error: failed to synthesize instance of type class
  Norm K.carrier
error: failed to synthesize instance of type class
  TopologicalSpace K.carrier
```

(`lean_check` で実測、0.12 秒。`import-audit --plan-all` はこの依存を検出できない
——コードに現れる**名前**しか追わないため。★これは見立てが外れた点として記録する。)

`ABC3.Found.PGC.LocalFieldNorm` が推移的に引く `Found` の本は
`ABC3.Found.ResidueFieldFinite` だけなので、この import で「引けなくなる」
`Found` の本は 200 本から数本に減る。
-/

namespace ABC3.Skeleton.PGC

open ABC3.Meta ABC3.Interface.PGC ABC3.Found.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ## Corollary 3.1 が語る対象 -/

/-- `RamificationFiltration`(`Interface/`)から `K` 1つ分の `FilteredGroup` を作る。

★`ABC3.Interface.PGC.RamificationFiltration.filt`(`Skeleton/PGC/Section4Defs.lean`)
および `ABC3.Found.PGC.filtOf`(`Found/PGC/RamificationNaturality.lean`)と
**同じ構造体リテラル**。import の向きの都合で 3 箇所に置いてあるが、
小さな構造の詰め替えなので複製のコストは小さい。 -/
noncomputable def filteredGroupOf (RF : RamificationFiltration p) (K : PAdicLocalField p) :
    FilteredGroup :=
  { G := K.absGal, Gv := RF.Gv K, isClosed := RF.isClosed K, isNormal := RF.isNormal K,
    antitone := RF.antitone K }

/-- 台帳の付随宣言(橋渡しの `def` であり、原典の項目そのものではない)。
`ABC3.Interface.PGC.RamificationFiltration.filt.src`(`Section4Defs.lean`)と同じ位置づけ。 -/
def filteredGroupOf.src : Source :=
  { paper := "pGC", pdfPage := 6, item := "Section 3 (filteredGroupOf)", sectionId := "cor-3-1" }

/-! ## Definition 3.2 -/

/-- **[pGC] Definition 3.2** — uniformizing な E[Γ_K]-加群。`sorry` 無し(純粋な定義)。

原文 (pGC p.6):
> We shall call the E[Γ_K]-module V uniformizing if the restriction of ρ_V to some open
> subgroup I of U_K (⊆ Γ^a_K^b) is the morphism I → E× induced by restricting some morphism
> of fields K → E to I ⊆ U_K ⊆ K.

## 形式化上の簡略化(逸脱として記録)

原文の「ρ_V は U_K ⊆ Γ_K^ab 上に制限できる」は、**局所類体論の相互律**
Γ_K^ab ≅ (K^×)^(§1 Proposition 1.2 の論拠、mathlib 不在)を暗黙に経由する——
U_K は本来 K^× の部分群であり、Γ_K の部分群ではない。

ここでは相互律を明示的なパラメータ `toGal : U_K → Γ_K`(未構築の辞書、まだ本物には
なっていない)として受け取ることで、定義全体を条件付きで well-typed にする
——`toGal` が本物になれば(Interface が実装されれば)この定義はそのまま使える。

## ★★★2026-09-05: 「開**部分群**」の条件が落ちていた(修理)

原文は "some open **subgroup** I of U_K" と言っているが、以前の形は
`I : Set K.carrier` に `IsOpen I` と `I ⊆ U_K` しか課しておらず、
**`I = ∅` が許されていた**。空集合では `∀ x ∈ I, …` が空虚に真になるので、
`E := K.carrier`・`ι := id` と取れば `IsUniformizing` は**どんな `ρ` でも
常に成り立つ**——定義が内容を失っていた
(`Check/PGC/Def32Degenerate.lean::isUniformizingOld_trivial`、`sorry` 無し)。

**修理**: 原文どおり `I` が `U_K` の**部分群**であること
(`1 ∈ I`・積で閉じる・逆元で閉じる)を課した。これで定義は識別力を持つ
——例えば自明な表現 `ρ = 1` は uniformizing に**ならない**
(`Check/PGC/Def32Degenerate.lean::not_isUniformizing_one`)。
`Lemma 4.1`(`Section4.lean`)が `1 ∈ I` を課していたのと同じ形。 -/
def IsUniformizing (K : PAdicLocalField p) (E : Type*) [Field E] [Algebra ℚ_[p] E]
    (toGal : {x : K.carrier // ‖x‖ = (1 : ℝ)} → K.absGal) (ρ : K.absGal →* Eˣ) : Prop :=
  ∃ (I : Set K.carrier) (hIU : I ⊆ {x : K.carrier | ‖x‖ = 1}) (_hopen : IsOpen I)
    (_hone : (1 : K.carrier) ∈ I)
    (_hmul : ∀ a ∈ I, ∀ b ∈ I, a * b ∈ I)
    (_hinv : ∀ a ∈ I, a⁻¹ ∈ I)
    (ι : K.carrier →+* E), ∀ x (hx : x ∈ I), ((ρ (toGal ⟨x, hIU hx⟩) : Eˣ) : E) = ι x

def IsUniformizing.src : Source :=
  { paper := "pGC", pdfPage := 6, item := "Definition 3.2", sectionId := "def-3-2" }

end ABC3.Skeleton.PGC
