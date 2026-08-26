/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.NormDL
import ABC3.Found.Divisor.SchemeCartier

/-!
# `K`-`Q`-Cartier と `Φ(L)^gp`(鎖 `normalize` の `KQC`、鎖 `sec6items` の `ex61-phigp`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109。

原文 (FrdI p.109):
> If for every Spec(L) ∈Ob(B(G)0) [cf. §0], every prime divisor of DL is Q-Cartier,

原文 (FrdI p.109):
> then we shall say that DK is K-Q-Cartier. In the following, we shall assume that

## ★★原文は仮定として置く

原文は「`D_K` が `K`-`Q`-Cartier である」を**仮定**として置く
(「we shall assume that `D_K` is `K`-`Q`-Cartier」)。★我々も**仮定として書く**。

★★その仮定が効くのは 1 か所だけである —— `Φ(L)^gp` にあたる部分群
(台が `D_L` に入る Cartier 因子)が **`IsQCartierSubgroup`** を満たすこと。
★これで `Found/Divisor/CartierMonoid.lean` 以降の単系論(perf-factorial まで)が
そのまま動く。

## ★本ファイルで閉じること

| 定義 | 中身 |
|---|---|
| `DLEmb` | `D_L ↪ V[L] の素因子` |
| `cartierOnDL` | ★**台が `D_L` に入る Cartier 因子の群**(= `Φ(L)^gp`) |
| `IsKQCartier` | ★原文の `K`-`Q`-Cartier(仮定) |
| `isQCartierSubgroup_cartierOnDL` | ★★**仮定 ⟹ `IsQCartierSubgroup`** |
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Found.FrdI

universe u

variable (V : Scheme.{u}) [IsIntegral V] [IsLocallyNoetherian V] {Kbar : Type u} [Field Kbar]
  [Algebra V.functionField Kbar]

/-! ## ★1. `D_L` から `V[L]` の素因子への埋め込み -/

/-- ★`D_L ↪ V[L] の素因子`。 -/
def DLEmb (DK : Set (PrimeDivisorPt V)) (L : FinSub V.functionField Kbar) :
    (DLSet V DK L) ↪ PrimeDivisorPt (normObj V L) :=
  Function.Embedding.subtype _

/-- ★台が `D_L` に入る Weil 因子を `V[L]` の Weil 因子として見る。 -/
noncomputable def toWeilOnDL (DK : Set (PrimeDivisorPt V)) (L : FinSub V.functionField Kbar) :
    ((DLSet V DK L) →₀ ℤ) →+ WeilDiv (normObj V L) where
  toFun d := Finsupp.embDomain (DLEmb V DK L) d
  map_zero' := by simp
  map_add' _ _ := by simp [Finsupp.embDomain_add]

/-! ## ★2. `Φ(L)^gp` —— 台が `D_L` に入る Cartier 因子 -/

/-- ★★**`Φ(L)^gp`** —— 台が `D_L` に入る Cartier 因子の群。

★`Found/Divisor/CartierMonoid.lean` の `Γ : AddSubgroup (S →₀ ℤ)` にそのまま入る形。 -/
noncomputable def cartierOnDL (DK : Set (PrimeDivisorPt V)) (L : FinSub V.functionField Kbar)
    [IsLocallyNoetherian (normObj V L)] (hnorm : IsNormalScheme (normObj V L)) :
    AddSubgroup ((DLSet V DK L) →₀ ℤ) :=
  (cartierSubgroup hnorm).comap (toWeilOnDL V DK L)

omit [IsLocallyNoetherian V] in
@[simp] theorem mem_cartierOnDL (DK : Set (PrimeDivisorPt V)) (L : FinSub V.functionField Kbar)
    [IsLocallyNoetherian (normObj V L)] (hnorm : IsNormalScheme (normObj V L))
    (d : (DLSet V DK L) →₀ ℤ) :
    d ∈ cartierOnDL V DK L hnorm ↔ IsCartierDiv hnorm (toWeilOnDL V DK L d) := Iff.rfl

/-! ## ★3. `K`-`Q`-Cartier -/

/-- ★★★**原文の `K`-`Q`-Cartier**(仮定)——
どの `L` についても `D_L` の各素因子が `Q`-Cartier。

原文 (FrdI p.109):
> then we shall say that DK is K-Q-Cartier. In the following, we shall assume that -/
def IsKQCartier (DK : Set (PrimeDivisorPt V))
    (hnorm : ∀ (L : FinSub V.functionField Kbar) [IsLocallyNoetherian (normObj V L)],
      IsNormalScheme (normObj V L)) : Prop :=
  ∀ (L : FinSub V.functionField Kbar) [IsLocallyNoetherian (normObj V L)]
    (s : DLSet V DK L),
    IsQCartierDiv (hnorm L) (Finsupp.single (s : PrimeDivisorPt (normObj V L)) 1)

omit [IsLocallyNoetherian V] in
/-- ★★★★**`K`-`Q`-Cartier ⟹ `Φ(L)^gp` は `IsQCartierSubgroup`**。

★★これで `CartierMonoid.lean` 以降の単系論(perf-factorial まで)がそのまま動く。 -/
theorem isQCartierSubgroup_cartierOnDL (DK : Set (PrimeDivisorPt V))
    (hnorm : ∀ (L : FinSub V.functionField Kbar) [IsLocallyNoetherian (normObj V L)],
      IsNormalScheme (normObj V L))
    (h : IsKQCartier V DK hnorm) (L : FinSub V.functionField Kbar)
    [IsLocallyNoetherian (normObj V L)] :
    IsQCartierSubgroup (cartierOnDL V DK L (hnorm L)) := by
  intro s
  obtain ⟨n, hn, hc⟩ := h L s
  refine ⟨n, hn, ?_⟩
  show IsCartierDiv (hnorm L) (toWeilOnDL V DK L (Finsupp.single s (n : ℤ)))
  have hemb : toWeilOnDL V DK L (Finsupp.single s (n : ℤ))
      = Finsupp.single (s : PrimeDivisorPt (normObj V L)) (n : ℤ) := by
    show Finsupp.embDomain (DLEmb V DK L) (Finsupp.single s (n : ℤ)) = _
    rw [Finsupp.embDomain_single]
    rfl
  rw [hemb]
  have hsm : (n : ℤ) • (Finsupp.single (s : PrimeDivisorPt (normObj V L)) (1 : ℤ))
      = Finsupp.single (s : PrimeDivisorPt (normObj V L)) (n : ℤ) := by
    rw [Finsupp.smul_single, smul_eq_mul, mul_one]
  have hnat : (n • (Finsupp.single (s : PrimeDivisorPt (normObj V L)) (1 : ℤ))
      : WeilDiv (normObj V L))
      = Finsupp.single (s : PrimeDivisorPt (normObj V L)) (n : ℤ) := by
    rw [← hsm, ← Nat.cast_smul_eq_nsmul ℤ]
  rw [← hnat]
  exact hc

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.1` の `K`-`Q`-Cartier と `Φ(L)^gp`。 -/
def IsKQCartier.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — K-Q-Cartier と Φ(L)^gp(台が D_L に入る Cartier 因子)",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
