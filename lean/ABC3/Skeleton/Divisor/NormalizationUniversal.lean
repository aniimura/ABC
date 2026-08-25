/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Skeleton.Divisor.Normalization
import ABC3.Found.Divisor.NormalSections

/-!
# 正規化の普遍性(正規スキームからの支配射)——`Theorem 6.2, (i)` の壁(`Skeleton`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.110。

原文 (FrdI p.110):
> in CV,K,DK may be thought of as consisting of the following data: (a) a morphism

## ★★なぜ要るか —— `thm62-i-pull` が要求するもの

`Theorem 6.2, (i)` は 2 つの幾何のデータ `(V₁, K̄₁, D_{K₁})`、`(V₂, K̄₂, D_{K₂})` の
間の射(仮定 (a)(b)(c))から `Ψ : C₁ → C₂` を作る。我々の側では

* 底の関手 `𝒟₁ → 𝒟₂` …… **在庫**(`thm62-i-Dfun`、`compFunctor`)
* `Φ₁ → Φ₂ ∘ ψ` …… **在庫**(`cartierPullback`)を `V₂ → V₁` に沿って当てる
* `B₁ → B₂ ∘ ψ` …… **在庫**(`normFFUnits`)

★★残るのは**因子を引き戻す土台**である。原文の `ψ : V₂ → V₁` は
**支配射であるだけ**(有限でも整でもない)なので、因子を引き戻すには

```
V₂[K₂(L)] ⟶ V₁[L]
```

が要る。★これは「**正規整スキームからの支配射は相対正規化を経由する**」であり、
★★**mathlib の `Scheme.Hom.normalizationDesc` では出ない**。

## ★★★★★測定 —— `normalizationDesc` は向きが逆である(2026-08-25 に再確認)

mathlib の普遍性は

```
normalizationDesc (f : X ⟶ Y) (f₁ : X ⟶ T) (f₂ : T ⟶ Y) [IsIntegralHom f₂]
    (H : f = f₁ ≫ f₂) : f.normalization ⟶ T
```

すなわち「**すでに `T` を経由する分解が与えられていれば**、正規化がそこへ降りる」
という形である。★我々が要るのは逆で、「`W` が正規で `W → Y` が支配的、
`L ⊆ K(W)` ならば `W → Y[L]` **を作る**」ほうである。
★★`f₁ : X ⟶ T` を先に要求するので、**そのままでは循環する**。

## ★中身(紙の上)

`Y[L] = Spec` of(`𝒪_Y` の `L` における整閉包)なので、`Y` 上の射 `W ⟶ Y[L]` は
`𝒪_Y`-代数の射 `IntegralClosure(𝒪_Y, L) → 𝒪_W` に対応する。
`W` が正規で `L ⊆ K(W)` なら、`𝒪_Y(U)` の `L` における整閉包は
`𝒪_W(ψ⁻¹U)` の中に入る(整な元は正規環の中で正則)。
★アフィン局所の主張はこれだけだが、mathlib の正規化は
**貼り合わせ**(`relativeGluingData`)で作られているので、
その貼り合わせをもう一度書くことになる。

★★★**(ii)(iii)(iv) は `V₁ = V₂` の場合**であり、そこでは `ψ` が Frobenius
(＝有限射)なので `normMap` と同じ筋で出る。**(i) の一般形だけ**がこの節点を要求する。
-/

namespace ABC3.Skeleton.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Meta

universe u

/-! ## ★1. 正規スキームからの支配射は相対正規化を経由する -/

open ABC3.Found.Divisor in
/-- ★★★★★★**正規化の普遍性(正規スキーム版)** ——
`W` が正規で `g : W ⟶ Y` が支配的、`L ⊆ K(W)` が `Y` の上で両立するなら、
`W` は `Y` の `L` における正規化 `f.normalization` を経由する。

## ★★★★測定の訂正(2026-08-25)—— 仮定の向きが逆だった

旧版は「`f` が `g` を経由する」(`e : SpecL ⟶ W`, `f = e ≫ g`)と書いていたが、
★**これは `K(W) ⊆ L` を意味していて、原文の `L ⊆ K(W)` の逆である**。
`Example 6.1` では `V[L]` の函数体が `L` で、`Theorem 6.2, (i)` が要求するのは
`V₂[K₂(L)] ⟶ V₁[L]`(函数体は `L₂ ⊇ L`)なので、正しい向きは

```
ε : Spec K(W) ⟶ SpecL   ただし  ε ≫ f = (生成点からの射) ≫ g
```

である。★仮定を直した。

## ★★★★★核は閉じた(2026-08-25)

数学の中身は「**正規スキームでは整な有理函数は切断**」1 本で、
`Found/Divisor/NormalSections.lean` の `mem_range_germ_of_forall_isIntegral` として
**`sorry` 無しで閉じた**。`Y[L]|_U = Spec(integralClosure Γ(Y,U) L)` の切断を
`Γ(W, g⁻¹U)` へ送る環準同型はこれで作れる。

★★残るのは**射の貼り合わせ**(`Scheme.OpenCover.glueMorphisms` で
`g⁻¹U ⟶ f.normalization` を貼る)という**配管**である。新しい数学は要らない。 -/
theorem exists_toNormalizationIn {W Y SpecL : Scheme.{u}} (f : SpecL ⟶ Y)
    [QuasiCompact f] [QuasiSeparated f] [IsIntegral Y] [IsIntegral W]
    (g : W ⟶ Y) [IsDominant g] (_hnorm : IsNormalScheme W)
    (_ε : Spec W.functionField ⟶ SpecL)
    (_hε : _ε ≫ f = W.fromSpecStalk (genericPoint (W : Type u)) ≫ g) :
    ∃ h : W ⟶ normalizationIn f, h ≫ f.fromNormalization = g := by
  sorry

/-- ★★**一意性** —— `Y` の上での射は高々 1 本(`fromNormalization` が mono だから)。

★`f.toNormalization` が支配的で `fromNormalization` が整なので、
`Y` 上の 2 本の射は一致する。 -/
theorem toNormalizationIn_unique {W Y SpecL : Scheme.{u}} (f : SpecL ⟶ Y)
    [QuasiCompact f] [QuasiSeparated f] [IsIntegral Y] [IsIntegral W]
    (g : W ⟶ Y) (h h' : W ⟶ normalizationIn f)
    (_hh : h ≫ f.fromNormalization = g) (_hh' : h' ≫ f.fromNormalization = g) :
    h = h' := by
  sorry

/-! ## ★2. 出口 —— 支配射に沿った `V[L]` の引き戻し -/

/-- ★★★★**`Theorem 6.2, (i)` が要求する射** ——
支配射 `ψ : V₂ ⟶ V₁` と体の埋め込み `L ↪ K₂(L)` から
`V₂[K₂(L)] ⟶ V₁[L]` を作る。

★★これがあれば `cartierPullback`(在庫)がそのまま当たり、
`thm62-i-pull` は**在庫の組み立て**になる(新しい数学は無い)。 -/
theorem exists_normalizationIn_map {V₁ V₂ SpecL SpecL₂ : Scheme.{u}}
    (f₁ : SpecL ⟶ V₁) (f₂ : SpecL₂ ⟶ V₂)
    [QuasiCompact f₁] [QuasiSeparated f₁] [QuasiCompact f₂] [QuasiSeparated f₂]
    [IsIntegral V₁] [IsIntegral V₂]
    (ψ : V₂ ⟶ V₁) [IsDominant ψ]
    (_hcompat : ∃ ε : SpecL₂ ⟶ SpecL, f₂ ≫ ψ = ε ≫ f₁) :
    ∃ Ψ : normalizationIn f₂ ⟶ normalizationIn f₁,
      Ψ ≫ f₁.fromNormalization = f₂.fromNormalization ≫ ψ := by
  sorry

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def exists_toNormalizationIn.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 正規スキームからの支配射は相対正規化を経由する",
    sectionId := "frdi-thm-6-2" }

/-- ★★★**`Theorem 6.2` の条なし `.src` を止めている 2 本のうちの 1 本**
(もう 1 本は `thm62-iv-slim`)。 -/
def exists_toNormalizationIn.needs : List ProofObligation :=
  [ .citation "[mathlib]" "Scheme.Hom.normalization / toNormalization / fromNormalization"
      (.inMathlib "AlgebraicGeometry.Scheme.Hom.normalization") 110,
    .citation "[mathlib]" "Scheme.Hom.normalizationDesc(★向きが逆。分解が先に要る)"
      (.inMathlib "AlgebraicGeometry.Scheme.Hom.normalizationDesc") 110,
    .citation "[mathlib]" "正規スキームからの支配射が正規化を経由すること"
      (.absent "Mathlib/AlgebraicGeometry/Normalization.lean の宣言を列挙(2026-08-25)。normalizationDesc / normalization.hom_ext はいずれも「分解が与えられていれば降りる」向きで、作る向きは無い") 110,
    .citation "[mathlib]" "IsIntegrallyClosed / integralClosure(アフィン局所の中身)"
      (.inMathlib "IsIntegrallyClosed") 110,
    .citation "[ABC3]" "★数学の中身(正規スキームでは整な有理函数は切断、sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.mem_range_germ_of_forall_isIntegral") 110,
    .citation "[ABC3]" "開集合の上での切断の貼り合わせ(sorry 無し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.exists_sectionOn") 110,
    .derivation
      "アフィン局所で「Γ(Y,U) の L における整閉包は Γ(W, ψ⁻¹U) に入る」——これは閉じた。残るのは Scheme.OpenCover.glueMorphisms による射の貼り合わせ" 110,
    .implicitStep
      "★原文は仮定 (a) から因子と有理函数の引き戻しを「(i)」の一言で置いている" 110 ]

def exists_normalizationIn_map.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 支配射に沿った V[L] の引き戻し",
    sectionId := "frdi-thm-6-2" }

def exists_normalizationIn_map.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_toNormalizationIn"
      (.inProject "ABC3" "ABC3.Skeleton.Divisor.exists_toNormalizationIn") 110,
    .citation "[ABC3]" "cartierPullback(支配射に沿った Cartier 因子の引き戻し)"
      (.inProject "ABC3" "ABC3.Found.Divisor.pullbackCartier") 110,
    .derivation "V₂[K₂(L)] は V₂ の上で正規なので、上の普遍性を V₁ の上で当てる" 110 ]

def toNormalizationIn_unique.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Theorem 6.2, (i) — 正規化への射は Y 上で一意",
    sectionId := "frdi-thm-6-2" }

def toNormalizationIn_unique.needs : List ProofObligation :=
  [ .citation "[mathlib]" "Scheme.Hom.normalization.hom_ext(正規化への 2 射の一致)"
      (.inMathlib "AlgebraicGeometry.Scheme.Hom.normalization.hom_ext") 110,
    .derivation "fromNormalization が整で toNormalization が支配的なので、Y 上の 2 射は一致する" 110 ]

end ABC3.Skeleton.Divisor
