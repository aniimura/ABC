/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop110

/-!
# sub-automorphism の不変量

★★[FrdI] `Proposition 1.6, (vi)` の `Aut^sub-ample` が
圏論的ファイバー積へ降りるかどうかを判定するために、
**sub-automorphism `α` の `degFr` と `Div` が何であるか**を先に測る。

`α ∈ End_𝒞(A)` が sub-automorphism であるとは、
ある `φ : B ⟶ A` と同型 `β ∈ Aut_𝒞(B)` が `β ≫ φ = φ ≫ α` を満たすことである(§0)。

## ★測った結果

1. **`degFr α = 1`**(`degFr_of_isSubAutomorphism`)。
   両辺の次数を取ると `degFr φ · 1 = degFr α · degFr φ` になり、`ℕ≥1` で消去できる。
2. **`Φ(Base β)(Div φ) = Φ(Base φ)(Div α) + Div φ`**
   (`div_relation_of_isSubAutomorphism`)。
   ★すなわち「`Φ(Base β)`(これは `Φ(Base B)` の**自己同型**)で `Div φ` を動かすと
   `Φ(Base φ)(Div α)` だけ増える」。
3. したがって `Φ(Base β)` が `Div φ` を動かさないなら
   **`Φ(Base φ)(Div α) = 0`**(`div_pullback_eq_zero_of_isSubAutomorphism`)。

## ★★`Div α = 0` は**言えない**(記録)

★2 の等式は `Div α = 0` を**強制しない**。`Φ(Base β)` が真に動かす例が作れる:

`𝒟` を対象 1 個・`End = ℤ`(すべて同型)の圏、`Φ(d) = ℝ≥0` で
`Φ.map n =`「`2ⁿ` 倍」、`B(d) = ℝ`(group-like)、`Div_B = id` とする。
`Theorem 5.2` の仮定(`Φ` divisorial・`B` group-like・`𝒟` connected かつ
totally epimorphic)はすべて満たすので model Frobenioid ができる。
そこで `Base β = b > 0`、`Base φ = ψ`、`Div φ = c > 0` と置き、
`Div α := (2^b − 1)·c / 2^ψ > 0` とすると 4 成分すべてが合う。

★**これが `Proposition 1.6, (vi)` の `Aut^sub-ample` を止めている理由である** ——
圏論的ファイバー積 `𝒞′ = 𝒞 ×_𝒟 𝒟′` で証人を作るには、
`𝒟′` の側の証人 `(e′, ψ′, β′)` に沿って `A` を引き戻した `χ : B ⟶ A` を取り、
引き戻しの普遍性で `γ₀`(`Base γ₀ = G β′`、`γ₀ ≫ χ = χ ≫ φ₀`)を得るが、
`Div γ₀ = Φ(Base χ)(Div φ₀)` であり、★これが `0` になる保証が無い。
`Div φ₀ = 0` なら済むが、上のとおり**それは言えない**。

★★上の 3 は「`Φ(Base β)` が `Div φ` を動かさない」ときの話で、
我々が要るのは `Φ(Base χ)(Div φ₀) = 0` —— **`χ` は証人の `φ` とは別の射**である。
したがって 3 では埋まらない。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-- ★★**sub-automorphism の Frobenius 次数は 1**。

★`β ≫ φ = φ ≫ α` の両辺の次数を取ると `degFr φ · 1 = degFr α · degFr φ` になり、
`degFr` は `ℕ≥1` の値なので消去できる。 -/
theorem degFr_of_isSubAutomorphism {A : C} {α : End A} (h : IsSubAutomorphism α) :
    P.degFr (α : A ⟶ A) = 1 := by
  obtain ⟨B, φ, β, hβ, heq⟩ := h
  have hd := congrArg P.degFr heq
  rw [P.degFr_comp, P.degFr_comp, degFr_of_isIso P (β : B ⟶ B), mul_one] at hd
  have hd2 : P.degFr φ * 1 = P.degFr φ * P.degFr (α : A ⟶ A) := by
    rw [mul_one, mul_comm]
    exact hd
  exact (mul_left_cancel hd2).symm

/-- ★★★**sub-automorphism の `Div` が満たす等式**。

`β ≫ φ = φ ≫ α` の両辺の `Div` を取ると、`β` が同型で `Div β = 0`、
`degFr α = 1` なので

  `Φ(Base β)(Div φ) = Φ(Base φ)(Div α) + Div φ`

になる。★`Φ(Base β)` は `Φ(Base B)` の**自己同型**(`β` が同型だから)なので、
これは「自己同型で `Div φ` を動かすと `Φ(Base φ)(Div α)` だけ増える」と言っている。 -/
theorem div_relation_of_isSubAutomorphism {A B : C} {α : End A} (φ : B ⟶ A) (β : End B)
    [hβ : IsIso (β : B ⟶ B)] (heq : (β : B ⟶ B) ≫ φ = φ ≫ (α : A ⟶ A)) :
    Φ.map (P.Base (β : B ⟶ B)) (P.Div φ)
      = Φ.map (P.Base φ) (P.Div (α : A ⟶ A)) + P.Div φ := by
  have hdb : P.Div (β : B ⟶ B) = 0 := isIsometric_of_isIso P (β : B ⟶ B)
  have hda : P.degFr (α : A ⟶ A) = 1 :=
    degFr_of_isSubAutomorphism P ⟨B, φ, β, hβ, heq⟩
  have h := congrArg P.Div heq
  rw [P.Div_comp, P.Div_comp, hdb, smul_zero, add_zero, hda] at h
  simpa using h

/-- ★★**証人の底自己同型が `Φ` に自明に作用するなら、sub-automorphism は
証人に沿って isometric** —— `Φ` が integral(消約的)なら消去できる。 -/
theorem div_pullback_eq_zero_of_isSubAutomorphism {A B : C} {α : End A} (φ : B ⟶ A)
    (β : End B) [hβ : IsIso (β : B ⟶ B)]
    (heq : (β : B ⟶ B) ≫ φ = φ ≫ (α : A ⟶ A))
    (htriv : Φ.map (P.Base (β : B ⟶ B)) (P.Div φ) = P.Div φ)
    (hcancel : ∀ x y z : Φ.val ((P.toElem.obj B).base), x + z = y + z → x = y) :
    Φ.map (P.Base φ) (P.Div (α : A ⟶ A)) = 0 := by
  have h := div_relation_of_isSubAutomorphism P φ β heq
  rw [htriv] at h
  exact hcancel _ _ (P.Div φ) (by rw [zero_add]; exact h.symm)

/-! ### ★出典の紐付け -/

/-- ★locator —— §0 の `sub-automorphism` の不変量。 -/
def degFr_of_isSubAutomorphism.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 14,
    item := "§0 sub-automorphism — Frobenius 次数は 1",
    sectionId := "frdi-s0-subaut" }

end ABC3.Found.FrdI
