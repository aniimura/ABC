import ABC3.Found.GaloisRep.MulByNTranscendental

/-!
# Galois (G5) 第 325 ブロック —— **★★★★★★`hfix` に要る 2 つの部品**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★道筋を測り直した——`[F(E):F(x)] = 2` は要らない

(G5) に残る `hfix`(`F(E)^{E[n]} = [n]^*F(E)`)は `deg[n] = n²` である。
★当初の筋は

    [F(E) : μF(E)] = [F(x) : F(x_n)]   ← `[F(E):F(x)] = 2` を μ の両側で使う

だったが、**`μ` に沿った次数の輸送**が重い。★★測り直すと、次の筋なら輸送が要らない:

    L := F(x_n, y_n)  と置く
    (a) x は L 上整で、最小多項式の次数は n² 以下   ← Φ_n − x_n·ΨSq_n がモニック(第 198)
    (b) y ∈ L(x)                                    ← 乗法公式の y 側の式
    (c) ゆえに L(x) = F(E)、つまり [F(E) : L] ≤ n²
    (d) L ⊆ Fix、[F(E) : Fix] = n²(第 196 の Artin)⟹ Fix = L
    (e) L ⊆ μF(E)(x_n = μx、y_n = μy)⟹ **hfix**

★★★★★**`[F(E):F(x)] = 2` も `μ` の輸送も要らない**——(b) が肩代わりする。

## ★★★★★★(b) の中身——乗法公式の `y` 側

第 52 の `MulOK` の第 3 成分

    (y_n − negY(x_n,y_n))·ΨSq_n(x)² = preΨ_{2n}(x)·(y − negY(x,y))

は `y − negY(x,y) = 2y + a₁x + a₃` なので、標数 0 なら `y` について解ける:

    y = ((2y_n + a₁x_n + a₃)·ΨSq_n(x)²/preΨ_{2n}(x) − a₁x − a₃)/2

★`preΨ_{2n}(x) ≠ 0` は `preΨ_{2n}` が非零多項式で `x` が超越的(第 116)だから。
★★右辺は `x_n, y_n, x` と `F` の元だけでできている——これが (b) である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `finrank_adjoin_le_of_monic_rel` | ★★★★★**(a)**——`z·B(x) = A(x)` から `x` の `L` 上の整性と次数評価 |
| `aeval_mem_intermediateField` | ★中間体は多項式の値で閉じている |
| `coordY_mem_of_mem` | ★★★★★★**(b)**——`x, x_n, y_n` を含む中間体は `y` も含む |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

/-! ## ★★★★★(a) 分子の次数が大きい関係からの整性と次数評価 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★**`z·B(x) = A(x)`、`A` モニックで `deg B < deg A` ならば
`x` は `z` を含む中間体 `L` 上整で、最小多項式の次数は `deg A` 以下**。

★`A − z·B` は `L` 係数のモニック多項式で次数は `deg A`(`deg B < deg A` だから)。 -/
theorem finrank_adjoin_le_of_monic_rel {F K : Type} [Field F] [Field K] [Algebra F K]
    (L : IntermediateField F K) (x z : K) (hz : z ∈ L)
    (A B : Polynomial F) (hA : A.Monic) (hdeg : B.natDegree < A.natDegree)
    (h : z * Polynomial.aeval x B = Polynomial.aeval x A) :
    IsIntegral L x ∧ (minpoly L x).natDegree ≤ A.natDegree := by
  set f : F →+* L := algebraMap F L with hf
  have hAm : (A.map f).Monic := hA.map f
  have hdm : (A.map f).natDegree = A.natDegree := hA.natDegree_map f
  have h1 : (Polynomial.C (⟨z, hz⟩ : L) * B.map f).natDegree ≤ B.natDegree := by
    refine Polynomial.natDegree_mul_le.trans ?_
    rw [Polynomial.natDegree_C, zero_add]
    exact Polynomial.natDegree_map_le
  have hlt : (Polynomial.C (⟨z, hz⟩ : L) * B.map f).natDegree < (A.map f).natDegree := by omega
  have hmonic : (A.map f - Polynomial.C (⟨z, hz⟩ : L) * B.map f).Monic := by
    rw [sub_eq_add_neg]
    refine Polynomial.Monic.add_of_left hAm ?_
    rw [Polynomial.degree_neg]
    exact Polynomial.degree_lt_degree hlt
  have hdeg2 : (A.map f - Polynomial.C (⟨z, hz⟩ : L) * B.map f).natDegree ≤ A.natDegree := by
    refine (Polynomial.natDegree_sub_le _ _).trans ?_
    simp only [max_le_iff]
    omega
  have heval : Polynomial.aeval x (A.map f - Polynomial.C (⟨z, hz⟩ : L) * B.map f) = 0 := by
    rw [map_sub, map_mul, Polynomial.aeval_C, Polynomial.aeval_map_algebraMap,
      Polynomial.aeval_map_algebraMap]
    show Polynomial.aeval x A - (z : K) * Polynomial.aeval x B = 0
    rw [h]; ring
  refine ⟨⟨_, hmonic, heval⟩, ?_⟩
  have hle := minpoly.degree_le_of_ne_zero L x hmonic.ne_zero heval
  have h2 := Polynomial.natDegree_le_natDegree hle
  omega

/-! ## ★★★★★★(b) `y` は `x, x_n, y_n` から出る -/

variable {F : Type} [Field F]

/-- ★中間体は多項式の値で閉じている。 -/
theorem aeval_mem_intermediateField {K : Type} [Field K] [Algebra F K]
    (M : IntermediateField F K) {t : K} (ht : t ∈ M) (p : Polynomial F) :
    Polynomial.eval₂ (algebraMap F K) t p ∈ M := by
  have h1 : Algebra.adjoin F ({t} : Set K) ≤ M.toSubalgebra := by
    rw [Algebra.adjoin_le_iff]
    rintro u hu
    rw [Set.mem_singleton_iff] at hu
    subst hu
    exact ht
  exact h1 (Polynomial.aeval_mem_adjoin_singleton F t)

set_option maxHeartbeats 1600000 in
/-- ★★★★★★**`x`・`x_n`・`y_n` を含む中間体は `y` も含む**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 52 の `MulOK` の `y` 側の式を `y` について解くだけである
(`preΨ_{2n}(x) ≠ 0` は `x` の超越性、割る `2` は標数 0)。 -/
theorem coordY_mem_of_mem [CharZero F] (W : WeierstrassCurve.Affine F)
    (M : IntermediateField F W.FunctionField) (n : ℕ) (hn : 1 ≤ n)
    {xn yn : W.FunctionField}
    (hy : (yn - (W.map (algebraMap F W.FunctionField)).negY xn yn)
            * Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W) (W.ΨSq (n : ℤ)) ^ 2
          = Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W) (W.preΨ (2 * (n : ℤ)))
            * (coordY W - (W.map (algebraMap F W.FunctionField)).negY (coordX W) (coordY W)))
    (hxn : xn ∈ M) (hyn : yn ∈ M) (hx : coordX W ∈ M) : coordY W ∈ M := by
  haveI : CharZero W.FunctionField :=
    charZero_of_injective_algebraMap (algebraMap F W.FunctionField).injective
  have hpre : Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W)
      (W.preΨ (2 * (n : ℤ))) ≠ 0 := by
    refine coordX_transcendental W (W.preΨ_ne_zero ?_)
    push_cast
    exact_mod_cast (by positivity : (0:ℚ) < 2 * (n:ℚ)) |>.ne'
  set a1 := algebraMap F W.FunctionField W.a₁ with ha1
  set a3 := algebraMap F W.FunctionField W.a₃ with ha3
  set S := Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W) (W.ΨSq (n : ℤ)) with hSdef
  set P := Polynomial.eval₂ (algebraMap F W.FunctionField) (coordX W)
    (W.preΨ (2 * (n : ℤ))) with hPdef
  have hneg1 : (W.map (algebraMap F W.FunctionField)).negY (coordX W) (coordY W)
      = -coordY W - a1 * coordX W - a3 := rfl
  have hneg2 : (W.map (algebraMap F W.FunctionField)).negY xn yn = -yn - a1 * xn - a3 := rfl
  rw [hneg1, hneg2] at hy
  have h2 : (2 : W.FunctionField) ≠ 0 := two_ne_zero
  have key : coordY W = ((2 * yn + a1 * xn + a3) * S ^ 2 / P - a1 * coordX W - a3) / 2 := by
    field_simp
    linear_combination -hy
  have hSm : S ∈ M := aeval_mem_intermediateField M hx _
  have hPm : P ∈ M := aeval_mem_intermediateField M hx _
  have h2m : (2 : W.FunctionField) ∈ M := by simp
  rw [key]
  refine M.div_mem (M.sub_mem (M.sub_mem (M.div_mem (M.mul_mem ?_ (M.pow_mem hSm 2)) hPm)
    (M.mul_mem ?_ hx)) ?_) h2m
  · exact M.add_mem (M.add_mem (M.mul_mem h2m hyn) (M.mul_mem (M.algebraMap_mem _) hxn))
      (M.algebraMap_mem _)
  · exact M.algebraMap_mem _
  · exact M.algebraMap_mem _

/-! ## ★出典の紐付け(`.src`) -/

def finrank_adjoin_le_of_monic_rel.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——deg[n] = n² の上からの評価)",
    sectionId := "genell-thm-3-8" }

def coordY_mem_of_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の非退化性——乗法公式の y 側)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
