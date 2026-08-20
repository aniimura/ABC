import ABC3.Found.GaloisRep.TwoTorsionCoord

/-!
# Galois (G1) 第 33 ブロック —— **★★★★★★`E[2]` は有限**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★G1 の第 2 欄(`torsion_finite`)の `n = 2` の場合

    {P | 2 • P = 0} は有限

★これが `torsion_finite` の**型そのもの**である——一般の `n` でも同じ骨で通る。

## ★★機構——3 つを繋ぐだけ

| 段 | 出所 |
|---|---|
| `2 • P = 0 ⟺ Ψ₂Sq(x_P) = 0` | ★第 31 |
| 2-捩れの `y` は `x` で決まる | ★第 32 |
| `Ψ₂Sq` の根は有限 | ★第 32(mathlib の `finite_setOf_isRoot`) |

★★★`x` 座標を取る写像 `xOf : Point → Option F` が**捩れ点の上で単射**であり、
像が `insert none (some '' 根の集合)` に入る——`Set.Finite.of_finite_image` で終わる。

## ★★一般の `n` へ

★同じ骨で `n • P = 0 ⟺ ΨSqₙ(x_P) = 0` と「`y` は高々 2 つ」を示せばよい。
★★`ΨSqₙ` が `X` だけの多項式なので根の有限性は同じ(mathlib の在庫)。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- 点の `x` 座標(無限遠点は `none`)。 -/
def xOf : W.toAffine.Point → Option F
  | .zero => none
  | .some x _ _ => some x

theorem finite_two_torsion (h2 : (2 : F) ≠ 0) :
    {P : W.toAffine.Point | (2 : ℕ) • P = 0}.Finite := by
  have hroot : ∀ P : W.toAffine.Point, (2 : ℕ) • P = 0 →
      xOf W P ∈ insert (none : Option F) ((fun x => (some x : Option F)) '' {x | W.Ψ₂Sq.IsRoot x}) := by
    rintro (_ | ⟨x, y, h⟩) hP
    · exact Set.mem_insert _ _
    · refine Set.mem_insert_of_mem _ ⟨x, ?_, rfl⟩
      exact (two_smul_eq_zero_iff W h).1 hP
  refine Set.Finite.of_finite_image (f := xOf W) ?_ ?_
  · have hfin : (insert (none : Option F)
        ((fun x => (some x : Option F)) '' {x | W.Ψ₂Sq.IsRoot x})).Finite :=
      Set.Finite.insert _ (Set.Finite.image _ (finite_psi2Sq_roots W h2))
    refine hfin.subset ?_
    rintro _ ⟨P, hP, rfl⟩
    exact hroot P hP
  · rintro (_ | ⟨x, y, h⟩) hP (_ | ⟨x', y', h'⟩) hP' hxx
    · rfl
    · exact absurd hxx (by simp [xOf])
    · exact absurd hxx (by simp [xOf])
    · simp only [xOf, Option.some.injEq] at hxx
      subst hxx
      have hy := (psi2Sq_eval_eq_zero_iff W h.left).1 ((two_smul_eq_zero_iff W h).1 hP)
      have hy' := (psi2Sq_eval_eq_zero_iff W h'.left).1 ((two_smul_eq_zero_iff W h').1 hP')
      have : y = y' := two_torsion_y_unique W h2 hy hy'
      subst this
      rfl


/-! ## ★出典の紐付け(`.src`) -/

def finite_two_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(E[2] の有限性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
