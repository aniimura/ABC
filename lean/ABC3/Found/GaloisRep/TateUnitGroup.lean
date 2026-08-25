import ABC3.Found.GaloisRep.TateMixedGroup

/-!
# Galois (G6) 第 273 ブロック —— **★★★★★★★★★単元 2 つの場合の群法則**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★到達点——場合分けが尽きた

基本領域 `0 ≤ v(u) < v(q)`(第 213)で正規化すると、`u·v·w = q` の 3 つ組は
`v(u) + v(v) + v(w) = v(q)` を満たす。★したがって `I` に入る個数は **1・2・3** のいずれか:

| `I` に入る個数 | 群法則 |
|---|---|
| 3 個(環帯だけ) | 第 261 |
| 2 個(単元 1 つ) | 第 272 |
| 1 個(単元 2 つ) | **本ブロック** |

★★★`v(u) = v(v) = 0` なら `v(w) = v(q)` なので `w` は環帯の**外縁**にある。
そこは `[w] = [q/(uv)] = [(uv)⁻¹]` であり、第 3 の点は `−P(uv)` である。
つまりこの場合が **`P(u) + P(v) = P(uv)`(単元どうしの準同型性)**にあたる。

## ★★★★★★`Λ` が `±` の分岐を潰す

`tate_inj_X_unit`(第 267)の第 2 の枝は `P(v) = −P(u)` である。これを潰すのに
**`Λ = Y/(X+Y)` を使う**:

    P(v) = −P(u)  ⟹  X' = X、Y' = −X−Y  ⟹  X'+Y' = −Y
                  ⟹  Λ(v) = (−X−Y)/(−Y) = 1/Λ(u)

★★★したがって `Λ(u)·Λ(v) = 1`。ところが `Λ(a) ≡ a mod I`(第 263)なので

    1 − u·v ∈ I   ⟹   1 − u·v は単元でない

★★★★**`IsUnit (1 − u·v)` が枝を殺す**——`hcp` にもとから入っている仮定である。
`u·v = 1` は `[v] = [u]⁻¹` すなわち第 3 の類が単位元、という退化にほかならない。

## ★★退化の 3 つ

相異性が壊れるのは 3 つの類のうち 2 つが一致するときで、それは

    u = v、  u·(u·v) = 1、  v·(u·v) = 1

に対応する(順に `[u] = [v]`、`[u] = [w]`、`[v] = [w]`)。
★本ブロックはこの 3 つを仮定で除く。**除いた場合(倍化)は補助母数で別に扱う**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `wOf_eq_of_mul` | ★相方は積で決まる |
| `isUnit_tateLambda` | ★★`Λ(a)` は単元 |
| `tateLambda_mul_eq_one_of_neg` | ★★★★★★`P(v) = −P(u)` なら `Λ(u)Λ(v) = 1` |
| `one_sub_mul_mem_of_neg` | ★★★★★★**`P(v) = −P(u)` なら `u·v ≡ 1`** |
| `tateXpair_ne_of_units` | ★★★★★★★★**単元どうしの座標の相異性** |
| `tate_points_add_eq_zero_two_units` | ★★★★★★★★★**単元 2 つの場合の群法則** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★下ごしらえ -/

theorem not_isUnit_of_mem [IsAdicComplete I R] [Nontrivial R] {x : R} (hx : x ∈ I) :
    ¬ IsUnit x := fun h => ne_of_isUnit_of_mem h hx rfl

theorem wOf_eq_of_mul {q a z : R} (ha : IsUnit a) (h : a * z = q) : wOf q a = z := by
  rw [wOf, ← h]
  calc a * z * Ring.inverse a = z * (a * Ring.inverse a) := by ring
    _ = z := by rw [Ring.mul_inverse_cancel _ ha, mul_one]

theorem ne_mul_self [Nontrivial R] {u v : R} (hu : IsUnit u) (hv1 : IsUnit (1 - v)) :
    u ≠ u * v := by
  intro h
  have hz : u * (1 - v) = 0 := by linear_combination h
  have h' := (hu.mul_right_eq_zero).1 hz
  rw [h'] at hv1
  exact not_isUnit_zero hv1

theorem isUnit_tateLambda [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (hw : w ∈ I)
    (ha : IsUnit a) (hu : IsUnit (1 - a)) : IsUnit (tateLambda a w q hq) := by
  have h := tateLambda_sub_mem a w q hq hw ha hu
  have he : tateLambda a w q hq = a + (tateLambda a w q hq - a) := by ring
  rw [he]
  exact isUnit_add_mem ha h

/-! ## ★★★★★★`Λ` が `±` の分岐を潰す -/

/-- ★★★★★★**`P(v) = −P(u)` なら `Λ(u)·Λ(v) = 1`**。

★`X' = X`・`Y' = −X−Y` なら `X'+Y' = −Y` なので `Λ(v) = (−X−Y)/(−Y) = 1/Λ(u)`。 -/
theorem tateLambda_mul_eq_one_of_neg [IsAdicComplete I R] (u v wu wv q : R) (hq : q ∈ I)
    (hwu : wu ∈ I) (hwv : wv ∈ I) (hu : IsUnit u) (hv : IsUnit v)
    (hu1 : IsUnit (1 - u)) (hv1 : IsUnit (1 - v))
    (hX : tateXpair v wv q hq = tateXpair u wu q hq)
    (hY : tateYpair v wv q hq = -tateXpair u wu q hq - tateYpair u wu q hq) :
    tateLambda u wu q hq * tateLambda v wv q hq = 1 := by
  set X := tateXpair u wu q hq with hXd
  set Y := tateYpair u wu q hq with hYd
  have hS : IsUnit (X + Y) := isUnit_tateXpair_add_tateYpair u wu q hq hwu hu hu1
  have hLu : tateLambda u wu q hq * (X + Y) = Y := by
    rw [tateLambda, mul_assoc, Ring.inverse_mul_cancel _ hS, mul_one]
  have hYu : IsUnit Y := by
    rw [← hLu]
    exact (isUnit_tateLambda u wu q hq hwu hu hu1).mul hS
  have hnY : IsUnit (-Y) := hYu.neg
  have hsum : tateXpair v wv q hq + tateYpair v wv q hq = -Y := by rw [hX, hY]; ring
  have hLv : tateLambda v wv q hq * (-Y) = -(X + Y) := by
    rw [tateLambda, hsum, mul_assoc, Ring.inverse_mul_cancel _ hnY, mul_one, hY]
    ring
  have hkey : ((X + Y) * (-Y)) * (tateLambda u wu q hq * tateLambda v wv q hq)
      = ((X + Y) * (-Y)) * 1 := by
    calc ((X + Y) * (-Y)) * (tateLambda u wu q hq * tateLambda v wv q hq)
        = (tateLambda u wu q hq * (X + Y)) * (tateLambda v wv q hq * (-Y)) := by ring
      _ = Y * (-(X + Y)) := by rw [hLu, hLv]
      _ = ((X + Y) * (-Y)) * 1 := by ring
  exact (hS.mul hnY).mul_left_cancel hkey

/-- ★★★★★★**`P(v) = −P(u)` なら `u·v ≡ 1`**——`Λ ≡ id` だから。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem one_sub_mul_mem_of_neg [IsAdicComplete I R] (u v wu wv q : R) (hq : q ∈ I)
    (hwu : wu ∈ I) (hwv : wv ∈ I) (hu : IsUnit u) (hv : IsUnit v)
    (hu1 : IsUnit (1 - u)) (hv1 : IsUnit (1 - v))
    (hX : tateXpair v wv q hq = tateXpair u wu q hq)
    (hY : tateYpair v wv q hq = -tateXpair u wu q hq - tateYpair u wu q hq) :
    1 - u * v ∈ I := by
  have h1 := tateLambda_mul_eq_one_of_neg u v wu wv q hq hwu hwv hu hv hu1 hv1 hX hY
  have hu' := tateLambda_sub_mem u wu q hq hwu hu hu1
  have hv' := tateLambda_sub_mem v wv q hq hwv hv hv1
  have hkey : 1 - u * v
      = u * (tateLambda v wv q hq - v) + (tateLambda u wu q hq - u) * tateLambda v wv q hq := by
    linear_combination -h1
  rw [hkey]
  exact Ideal.add_mem _ (Ideal.mul_mem_left _ _ hv') (Ideal.mul_mem_right _ _ hu')

/-- ★★★★★★★★**単元どうしの座標の相異性**——`u ≠ v` と `1 − u·v` が単元。 -/
theorem tateXpair_ne_of_units [IsAdicComplete I R] [IsDomain R] (u v wu wv q : R) (hq : q ∈ I)
    (hu : IsUnit u) (hv : IsUnit v) (hu1 : IsUnit (1 - u)) (hv1 : IsUnit (1 - v))
    (hwu : u * wu = q) (hwv : v * wv = q) (hwuI : wu ∈ I) (hwvI : wv ∈ I)
    (huv : u ≠ v) (huv1 : IsUnit (1 - u * v)) :
    tateXpair u wu q hq ≠ tateXpair v wv q hq := by
  intro heq
  have e1 : wOf q u = wu := wOf_eq_of_mul hu hwu
  have e2 : wOf q v = wv := wOf_eq_of_mul hv hwv
  rcases tate_inj_X_unit u v q hq hu hv hu1 hv1 (by rw [e1, e2]; exact heq) with h | h
  · exact huv h
  · rw [e1, e2] at h
    exact not_isUnit_of_mem
      (one_sub_mul_mem_of_neg u v wu wv q hq hwuI hwvI hu hv hu1 hv1 heq.symm h) huv1

/-! ## ★★★★★★★★★単元 2 つの場合の群法則 -/

theorem collUnits_two_units [IsAdicComplete I R] {u v w : R} (hu1 : IsUnit (1 - u))
    (hv1 : IsUnit (1 - v)) (hw : w ∈ I) (huv1 : IsUnit (1 - u * v)) :
    ∀ i, IsUnit (1 - collPts u v w i) := by
  intro i
  fin_cases i
  · exact hu1
  · exact hv1
  · exact isUnit_one_sub hw
  · exact isUnit_one_sub (Ideal.mul_mem_left _ _ hw)
  · exact isUnit_one_sub (Ideal.mul_mem_left _ _ hw)
  · exact huv1

section Group

variable {K : Type} [IsDomain R] [IsAdicComplete I R] [Field K] [DecidableEq K] [Algebra R K]

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★**単元 2 つの場合の群法則**——`P(u) + P(v) = P(u·v)` にあたる。

★退化の 3 つ(`u = v`、`u·(u·v) = 1`、`v·(u·v) = 1`)は仮定で除く。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_points_add_eq_zero_two_units (hinj : Function.Injective (algebraMap R K))
    (u v w q : R) (hq : q ∈ I) (hu : IsUnit u) (hv : IsUnit v)
    (hu1 : IsUnit (1 - u)) (hv1 : IsUnit (1 - v)) (hw : w ∈ I)
    (huvw : u * v * w = q) (huv : u ≠ v)
    (huv1 : IsUnit (1 - u * v)) (huuv : IsUnit (1 - u * (u * v)))
    (hvuv : IsUnit (1 - v * (u * v)))
    (n₁ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (algebraMap R K (tateXpair u (v * w) q hq)) (algebraMap R K (tateYpair u (v * w) q hq)))
    (n₂ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (algebraMap R K (tateXpair v (u * w) q hq)) (algebraMap R K (tateYpair v (u * w) q hq)))
    (n₃ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (algebraMap R K (tateXpair w (u * v) q hq)) (algebraMap R K (tateYpair w (u * v) q hq))) :
    Point.some _ _ n₁ + Point.some _ _ n₂ + Point.some _ _ n₃ = 0 := by
  have hvwI : v * w ∈ I := Ideal.mul_mem_left _ _ hw
  have huwI : u * w ∈ I := Ideal.mul_mem_left _ _ hw
  have huvU : IsUnit (u * v) := hu.mul hv
  have e1 : u * (v * w) = q := by rw [← huvw]; ring
  have e2 : v * (u * w) = q := by rw [← huvw]; ring
  have e3 : (u * v) * w = q := huvw
  have hsym : tateXpair w (u * v) q hq = tateXpair (u * v) w q hq := tateXpair_symm (u * v) w q hq
  have h12 := tateXpair_ne_of_units u v (v * w) (u * w) q hq hu hv hu1 hv1 e1 e2 hvwI huwI
    huv huv1
  have h13 := tateXpair_ne_of_units u (u * v) (v * w) w q hq hu huvU hu1 huv1 e1 e3 hvwI hw
    (ne_mul_self hu hv1) huuv
  have h23 := tateXpair_ne_of_units v (u * v) (u * w) w q hq hv huvU hv1 huv1 e2 e3 huwI hw
    (by rw [show u * v = v * u from mul_comm u v]; exact ne_mul_self hv hu1) hvuv
  rw [← hsym] at h13 h23
  have hd0 := TateCollinearSection.tate_collinear u v w q hq huvw
    (collUnits_two_units hu1 hv1 hw huv1)
  have hd := congrArg (algebraMap R K) hd0
  simp only [map_add, map_mul, map_sub, map_zero] at hd
  exact add_add_eq_zero_of_collDet n₁ n₂ n₃ hd (fun h => h12 (hinj h))
    (fun h => h13 (hinj h)) (fun h => h23 (hinj h))

end Group

/-! ## ★出典の紐付け(`.src`) -/

def one_sub_mul_mem_of_neg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——P(v) = -P(u) なら uv ≡ 1)",
    sectionId := "genell-def-3-3" }

def tate_points_add_eq_zero_two_units.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——単元 2 つの場合の群法則)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
