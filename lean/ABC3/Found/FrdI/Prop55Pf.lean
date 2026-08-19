/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def31Pf
import ABC3.Found.FrdI.Prop25iii

/-!
# [FrdI] Proposition 5.5, (i) の中身 —— `𝒞^pf` の遷移写像は `α ↦ α^m`

原文 (FrdI p.104):
> (i) If A ∈Ob(Cistr) maps to an object Apf ∈Ob(Cpf), then the natural functor

★原文 `(i)` は `𝒪^▷(A)^pf ≅ 𝒪^▷(A^pf)` を主張し、証明を
「Frobenius-trivial な `A` については、base-identity な Frobenius 型自己射を考え、
`𝒞` が Frobenius-normalized 型であることを使えば **immediately**」で畳んでいる。

★★**その "immediately" の中身がこれ**である。

## ★なぜこれで済むのか

`𝒞^pf` の射は `Definition 3.1, (iii)` の**余極限**で定義され、
その遷移写像は `frobTransport`(`Proposition 1.10, (i)` の存在と**一意性**)で与えられる。
`A` が Frobenius-trivial なら、添字は base-identity な Frobenius 型自己射 `ζ_m`
(次数 `m`)で走らせられる。★そこで `𝒪^▷(A)` の元 `α` に対し遷移写像を計算すると、
**Frobenius-normalized**(`ζ ≫ α^m = α ≫ ζ`)と一意性から

  `frobTransport ζ ζ α = α ^ m`

★★したがって `𝒪^▷(A^pf)` は「`𝒪^▷(A)` を `α ↦ α^m` で繋いだ余極限」、
すなわち **perfection `𝒪^▷(A)^pf`** になる。

## ★★残り(記録)

余極限そのものの同定(`ℕ≥1` で添字づけた族が `IdxPf A A` の中で **cofinal** であること、
および `Pf (Additive (𝒪^▷(A)))` との突き合わせ)は別の葉である
—— 依存グラフの鎖 `prop55` の `p55i`。
★本ファイルが与えるのは**その計算の核**である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P}

/-- ★★★★**`𝒞^pf` の遷移写像は `𝒪^▷(A)` の上で `α ↦ α^m`**。

★証明は `frobTransport` の**一意性**(`Proposition 1.10, (i)`)1 本:
Frobenius-normalized が `ζ ≫ α^m = α ≫ ζ` を与えるので、`α^m` が遷移写像に一致する。 -/
theorem frobTransport_otri_pow {A : C} (hfn : IsFrobeniusNormalized P A)
    (ζ : End A) (hζ : IsFrobeniusType P ((ζ : A ⟶ A))) (hζb : IsBaseIdentity P ζ)
    (α : End A) (hα : α ∈ OTri P A) :
    frobTransport (F := F) ((ζ : A ⟶ A)) hζ ((ζ : A ⟶ A)) hζ rfl ((α : A ⟶ A))
      = ((α ^ ((P.degFr ((ζ : A ⟶ A)) : ℕ+) : ℕ) : End A) : A ⟶ A) :=
  frobTransport_eq _ hζ _ hζ rfl _ _ (hfn ζ hζb α hα).symm

/-- ★遷移した先も `𝒪^▷(A)` に入る(部分単系だから)。 -/
theorem frobTransport_otri_mem {A : C} (hfn : IsFrobeniusNormalized P A)
    (ζ : End A) (hζ : IsFrobeniusType P ((ζ : A ⟶ A))) (hζb : IsBaseIdentity P ζ)
    (α : End A) (hα : α ∈ OTri P A) :
    (frobTransport (F := F) ((ζ : A ⟶ A)) hζ ((ζ : A ⟶ A)) hζ rfl ((α : A ⟶ A)) : End A)
      ∈ OTri P A := by
  rw [frobTransport_otri_pow hfn ζ hζ hζb α hα]
  exact (OTri P A).pow_mem hα _

/-- ★`𝒪^▷(A)` の可換性を `End A` の言葉で(`Commute` 版)。 -/
theorem otri_comm_end' {A : C} (hfn : IsFrobeniusNormalized P A) {x y : End A}
    (hx : x ∈ OTri P A) (hy : y ∈ OTri P A) : Commute x y :=
  congrArg Subtype.val (otri_mul_comm P hfn ⟨x, hx⟩ ⟨y, hy⟩)

/-- ★★**遷移写像は `𝒪^▷(A)` の単系準同型** —— `α ↦ α^m`。

★`𝒪^▷(A)` は Frobenius-normalized なら**可換**(`Proposition 2.5, (iii)`)なので、
`α ↦ α^m` は単系準同型になる。 -/
noncomputable def otriPowHom {A : C} (hfn : IsFrobeniusNormalized P A) (m : ℕ) :
    OTri P A →* OTri P A where
  toFun α := ⟨((α : End A) ^ m : End A), (OTri P A).pow_mem α.2 m⟩
  map_one' := Subtype.ext (one_pow m)
  map_mul' x y := Subtype.ext (by
    show ((x : End A) * (y : End A)) ^ m = ((x : End A) ^ m) * ((y : End A) ^ m)
    exact (otri_comm_end' hfn x.2 y.2).mul_pow m)

/-! ### ★出典の紐付け -/

/-- ★locator —— `Proposition 5.5, (i)` の「immediately」の中身
(★**条つき**: 余極限の同定そのものは未実装)。 -/
def frobTransport_otri_pow.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (i) — 𝒞^pf の遷移写像は 𝒪^▷(A) の上で α ↦ α^m",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
