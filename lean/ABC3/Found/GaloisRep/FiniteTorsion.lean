import ABC3.Found.GaloisRep.FiniteY

/-!
# Galois (G1) 第 35 ブロック —— **★★★★★★`torsion_finite` は乗法公式の片側だけに帰着する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★残る穴を 1 点に絞った

G1 の第 2 欄 `torsion_finite` は、次の**片側の含意**さえあれば出る:

    m • P = 0  ⟹  ΨSqₙ(x_P) = 0

★★逆向き(`ΨSqₙ(x) = 0 ⟹ m • P = 0`)は**要らない**——有限性には上界だけで足りる。

## ★★機構

| 段 | 出所 |
|---|---|
| `ΨSqₙ` の根は有限 | ★mathlib `ΨSq_ne_zero` + `finite_setOf_isRoot` |
| 1 つの `x` の上の `y` は有限 | ★第 34 |
| 座標を取る写像は単射 | ★`Point.some` の構成子 |

★★★`{(x,y) | ΨSqₙ(x) = 0 ∧ Equation x y}` が有限——根ごとに `y` を集めるだけ。

## ★★これで G1 の構図が確定した

| 欄 | 残り |
|---|---|
| `torsion_finite` | ★**乗法公式の片側のみ** |
| `structure_eq` | ★★乗法公式 + 分離性 + 次数(次数は mathlib の在庫) |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- 点の座標(無限遠点は `none`)。 -/
def xyOf : W.toAffine.Point → Option (F × F)
  | .zero => none
  | .some x y _ => some (x, y)

theorem finite_torsion_of_formula {m : ℕ} {n : ℤ} (hn : (n : F) ≠ 0)
    (hform : ∀ (x y : F) (h : W.toAffine.Nonsingular x y),
      m • (Point.some x y h) = 0 → (W.ΨSq n).IsRoot x) :
    {P : W.toAffine.Point | m • P = 0}.Finite := by
  classical
  have hroots : {x : F | (W.ΨSq n).IsRoot x}.Finite :=
    Polynomial.finite_setOf_isRoot (W.ΨSq_ne_zero hn)
  have hpairs : {xy : F × F | (W.ΨSq n).IsRoot xy.1 ∧ W.toAffine.Equation xy.1 xy.2}.Finite := by
    refine Set.Finite.subset (Set.Finite.biUnion hroots
      (fun x _ => Set.Finite.image (fun y => (x, y)) (finite_y_of_x W x))) ?_
    rintro ⟨x, y⟩ ⟨hx, hy⟩
    exact Set.mem_biUnion hx ⟨y, hy, rfl⟩
  refine Set.Finite.of_finite_image (f := xyOf W) ?_ ?_
  · refine Set.Finite.subset (Set.Finite.insert none
      (Set.Finite.image (fun xy => (some xy : Option (F × F))) hpairs)) ?_
    rintro _ ⟨P, hP, rfl⟩
    rcases P with _ | ⟨x, y, h⟩
    · exact Set.mem_insert _ _
    · exact Set.mem_insert_of_mem _ ⟨(x, y), ⟨hform x y h hP, h.left⟩, rfl⟩
  · rintro (_ | ⟨x, y, h⟩) hP (_ | ⟨x', y', h'⟩) hP' hxx
    · rfl
    · exact absurd hxx (by simp [xyOf])
    · exact absurd hxx (by simp [xyOf])
    · simp only [xyOf, Option.some.injEq, Prod.mk.injEq] at hxx
      obtain ⟨rfl, rfl⟩ := hxx
      rfl


/-! ## ★出典の紐付け(`.src`) -/

def finite_torsion_of_formula.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(torsion_finite は乗法公式の片側に帰着すること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
