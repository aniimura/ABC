import ABC3.Found.GaloisRep.DegBoundParts

/-!
# Galois (G5) 第 326 ブロック —— **★★★★★★`hfix` の帳簿の道具**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★第 325 の (c) と (d) を道具として切り出した

第 325 で `hfix` の道筋を

    L := F(x_n, y_n)
    (a) x は L 上整、最小多項式の次数 ≤ n²      ✅ 第 325
    (b) y ∈ L(x)                                 ✅ 第 325
    (c) ゆえに L(x) = F(E)                       ★本ブロック
    (d) L ⊆ Fix、[F(E) : Fix] = n² と挟む        ★本ブロック
    (e) L ⊆ μF(E) ⟹ hfix

と定めた。★本ブロックは (c) と (d) を、**曲線から独立した形**で用意する。

## ★★★★★(c) —— `x` と `y` を含む中間体は `⊤`

`F[W]` は `AdjoinRoot` なので `CoordinateRing.mk` が全射であり、`F[X][Y]` の元は
単項式の和である。★単項式 `C a · Y^k` の像は `a(x)·y^k` で、
中間体は多項式の値と積で閉じているから `M` に入る(`algebraMap_mk_C`)。
★★`F(E)` は `F[W]` の分数体だから、`IsFractionRing.div_surjective` で `⊤` になる。

## ★★★★★(d) —— 次数で挟んで中間体を一致させる

`L ≤ M`、`[K : L] ≤ N`、`[K : M] = N`、`0 < [K : L]` なら `L = M`。
★`IntermediateField.extendScalars` で `M` を `L` 上の中間体と見て、
`finrank L M · finrank M K = finrank L K` の**塔**に当てるだけである。
★★`finrank L M = 1` から `M = ⊥`、すなわち `M ⊆ L` が出る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `algebraMap_mk_C` | ★`mk (C p)` の像は `p(x)` |
| `intermediateField_eq_top_of_coord` | ★★★★★★**(c)** `x, y` を含む中間体は `⊤` |
| `finrank_le_of_adjoin_top` | ★★★★`L(x) = K` から `[K:L] ≤ deg(minpoly)` |
| `eq_of_finrank_bound` | ★★★★★**(d)** 次数で挟んで中間体を一致させる |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IntermediateField

variable {F : Type} [Field F]

/-! ## ★★★★★★(c) `x` と `y` を含む中間体は `⊤` -/

/-- ★`mk (C p)` の関数体での像は `p(x)` である。 -/
theorem algebraMap_mk_C (W : WeierstrassCurve.Affine F) (p : Polynomial F) :
    algebraMap W.CoordinateRing W.FunctionField (CoordinateRing.mk W (Polynomial.C p))
      = Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W) p := by
  have h : ((algebraMap W.CoordinateRing W.FunctionField).comp
      ((CoordinateRing.mk W).comp (Polynomial.C : Polynomial F →+* Polynomial (Polynomial F))))
      = Polynomial.eval₂RingHom (algebraMap F W.FunctionField) (coordX W) := by
    refine Polynomial.ringHom_ext ?_ ?_
    · intro a
      simp [coordX, genX]
      rfl
    · simp [coordX, genX]
  exact DFunLike.congr_fun h p

set_option maxHeartbeats 1600000 in
/-- ★★★★★★**`x` と `y` を含む中間体は `F(E)` 全体である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`F[W]` は `genX, genY` で生成され、`F(E)` はその分数体だから。 -/
theorem intermediateField_eq_top_of_coord (W : WeierstrassCurve.Affine F)
    (M : IntermediateField F W.FunctionField)
    (hx : coordX W ∈ M) (hy : coordY W ∈ M) : M = ⊤ := by
  have hmem : ∀ a : W.CoordinateRing,
      algebraMap W.CoordinateRing W.FunctionField a ∈ M := by
    intro a
    obtain ⟨q, rfl⟩ := Quotient.exists_rep a
    show algebraMap W.CoordinateRing W.FunctionField (CoordinateRing.mk W q) ∈ M
    induction q using Polynomial.induction_on' with
    | add p q hp hq =>
      rw [map_add, map_add]
      exact M.add_mem hp hq
    | monomial k a =>
      rw [← Polynomial.C_mul_X_pow_eq_monomial, map_mul, map_mul, map_pow, map_pow,
        algebraMap_mk_C]
      refine M.mul_mem (aeval_mem_intermediateField M hx _) (M.pow_mem ?_ k)
      show algebraMap W.CoordinateRing W.FunctionField (genY W) ∈ M
      exact hy
  refine eq_top_iff.2 (fun z _ => ?_)
  obtain ⟨a, b, hb, hab⟩ := IsFractionRing.div_surjective (A := W.CoordinateRing) z
  rw [← hab]
  exact M.div_mem (hmem a) (hmem b)

/-! ## ★★★★★(d) 次数で挟んで中間体を一致させる -/

set_option maxHeartbeats 1600000 in
/-- ★★★★`L(x) = K` なら `[K : L]` は最小多項式の次数に等しい。 -/
theorem finrank_le_of_adjoin_top {F K : Type} [Field F] [Field K] [Algebra F K]
    (L : IntermediateField F K) (x : K) (hint : IsIntegral L x) (N : ℕ)
    (hb : (minpoly L x).natDegree ≤ N)
    (htop : IntermediateField.restrictScalars F (L⟮x⟯) = ⊤) :
    Module.finrank L K ≤ N ∧ 0 < Module.finrank L K := by
  have h1 : (L⟮x⟯ : IntermediateField L K) = ⊤ := by
    refine IntermediateField.restrictScalars_injective F ?_
    rw [htop, IntermediateField.restrictScalars_top]
  have h2 : Module.finrank L ↥(L⟮x⟯ : IntermediateField L K) = (minpoly L x).natDegree :=
    IntermediateField.adjoin.finrank hint
  rw [h1, IntermediateField.finrank_top'] at h2
  have h3 : 0 < (minpoly L x).natDegree := minpoly.natDegree_pos hint
  omega

set_option maxHeartbeats 1600000 in
/-- ★★★★★**次数で挟んで中間体を一致させる**。

`L ≤ M`、`[K : L] ≤ N`、`[K : M] = N`、`0 < [K : L]` ならば `L = M`。

★`extendScalars` で `M` を `L` 上の中間体と見て、塔 `[M:L]·[K:M] = [K:L]` に当てる。 -/
theorem eq_of_finrank_bound {F K : Type} [Field F] [Field K] [Algebra F K]
    (L M : IntermediateField F K) (hLM : L ≤ M) (N : ℕ)
    (hpos : 0 < Module.finrank L K)
    (hL : Module.finrank L K ≤ N) (hM : Module.finrank M K = N) :
    L = M := by
  set M' : IntermediateField L K := IntermediateField.extendScalars hLM with hM'
  have htower : Module.finrank L M' * Module.finrank M' K = Module.finrank L K :=
    Module.finrank_mul_finrank L M' K
  have hMK : Module.finrank M' K = N := hM
  rw [hMK] at htower
  have h1 : Module.finrank L M' = 1 := by
    rcases Nat.eq_zero_or_pos (Module.finrank L M') with h | h
    · rw [h] at htower; omega
    · nlinarith [htower, hL, hpos, h]
  have hbot : M' = ⊥ := IntermediateField.finrank_eq_one_iff.1 h1
  refine le_antisymm hLM (fun m hm => ?_)
  have hmM' : m ∈ M' := by
    rw [hM']
    exact (IntermediateField.mem_extendScalars hLM).2 hm
  rw [hbot, IntermediateField.mem_bot] at hmM'
  obtain ⟨u, hu⟩ := hmM'
  rw [← hu]
  exact u.2

/-! ## ★出典の紐付け(`.src`) -/

def intermediateField_eq_top_of_coord.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——関数体は座標関数で生成される)",
    sectionId := "genell-thm-3-8" }

def eq_of_finrank_bound.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——次数で挟む)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
