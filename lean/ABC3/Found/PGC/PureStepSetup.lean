import ABC3.Found.PGC.GainedJumpFree
import Mathlib.FieldTheory.Galois.Basic
import Mathlib.Analysis.Normed.Unbundled.SpectralNorm
import ABC3.Found.PGC.LocalFieldNorm

/-!
# [pGC] ★★★★★★配られた配管 10 本 —— 4 本は落ちた、そして **`_of_conj_pure` は `k ≥ 1` で空虚**だった

持ち場は「`GainedJumpFree.exists_norm_sub_algebraMap_le_prod_axDecay_of_conj_pure` の
配管 10 本(`he hnormp heM hval hdeg htop hsep hsp hconj hne`)を落とす」であった。

## ★★★①配られた字面の検算 —— **`hne` は `1 ≤ k` で偽**である

`_of_conj_pure` は

```
τ : M →ₐ[F] M,  (minpoly F π).natDegree = p,  Algebra.adjoin F {π} = ⊤,
hne : (τ ^ p ^ k) π ≠ π
```

を仮定する。★後ろ 2 本は `[M : F] = p` を意味する。すると

* `M` は `F` 上有限次なので `τ` は全単射(`AlgHom.bijective`)、すなわち `τ ∈ Aut(M/F)`、
* Artin(`IntermediateField.finrank_fixedField_eq_card`)と塔の公式で
  `orderOf τ ∣ [M : F] = p`、
* ゆえに `τ ^ p = 1`、したがって `1 ≤ k` なら `τ ^ (p ^ k) = 1` で `π` を**動かせない**。

★★これを形式化したのが `algHom_pow_pow_apply_eq_self` / `conj_pure_absurd_of_one_le` である。
★★**`_of_conj_pure`(および `_of_conj_free` / `_of_jumps_free`)は `1 ≤ k` では
仮説が両立せず、空虚に真**である。使えるのは `k = 0` だけで、そのとき定数は
`∏_{j ∈ Icc 1 1} axDecay p j = axDecay p 1 = p^{1/(p−1)}` になる。

### ★なぜこうなったか(`GainedJumpSeq` との差)

`GainedJumpSeq.exists_norm_sub_algebraMap_le_prod_axDecay_lt` では
★`τ : M →+ M`(**加法準同型だけ**)であり、`F`-代数準同型なのは `σ = τ^{p^k}` の方だった。
これは正しい設定である —— 塔 `F_0 ⊂ F_1 ⊂ … ⊂ F_{k+1} = M` で `τ` は `Gal(M/F_0)` の
生成元、`F = F_k`、`τ^{p^k}` がちょうど `Gal(M/F)` を生成する。

★`GainedJumpFree` は `hstep` を仮説から**導出**するために `τ` を `M →ₐ[F] M` に**強めた**
(各層 `j < k` で `norm_algHom_sub_le_zpow_mul` を使うため `τ^{p^j}` が `F`-線形である必要がある)。
★★しかし層 `j` で要るのは「`F_j`-線形」であって「`F_k`-線形」ではない。
★この 1 段の取り違えで `k ≥ 1` が潰れている。

★★★**したがって「跳びの列も Hasse–Arf も要らなくなった」という結論は、
`k = 0`(= 巡回 `p` 次 1 段)についてのみ正しい。** `k = 0` の主張は
`MinpolyOrbitSplit.exists_norm_sub_algebraMap_le_axDecay_of_orbit` として**既に木に在る**。

## ★★②何を作ったか —— 抽象核と `WildStep`

| 節 | 内容 | 分岐・付値の語彙 |
|---|---|---|
| §1 | `[M : F] = n` ⟹ `τ ^ n = 1`(★Artin + 塔の公式だけ) | ★1 語も出ない |
| §2 | `_of_conj_pure` の仮説は `1 ≤ k` で矛盾 | ★1 語も出ない |
| §3 | スペクトルノルムなら `k`-代数準同型は等長(`hiso` の供給) | ノルムのみ |
| §4 | `structure WildStep` —— 設定を 1 つにまとめ、4 本を定理として出す | 具体層 |
| §5 | `PAdicLocalField` からの供給(`hnormp` / `hiso`) | 具体層 |

## ★★③配管 10 本の帰趨(★1 本ずつ)

| 仮説 | 帰趨 |
|---|---|
| `hsep`(分離性) | ★★**落ちた** —— `WildStep.hsep`(軌道が相異なる `p` 点) |
| `hsp`(分解) | ★★**落ちた** —— `WildStep.hsp` |
| `hconj`(全共役の距離) | ★★**落ちた** —— `WildStep.hconj` |
| `hne`(`τ^{p^k}π ≠ π`) | ★★**落ちた** —— `WildStep.hne`(`hjump` と `‖π‖ > 0` から) |
| `hnormp`(`‖p‖ = 1/p`) | ★★**落ちた**(p 進局所体では) —— `norm_natCast_p_eq_inv` |
| `hdeg` / `htop` | ★残る。**`M = F(π)` が `p` 次**という設定そのもの |
| `heM` / `he` | ★残る。**`v_M(p) = p·e`**(絶対分岐指数)という設定そのもの |
| `hval` | ★残る。**`M/F` が全分岐**というノルム言語での中身 |

さらに `GainedJumpFree` にはなかったが `MinpolyOrbitSplit` 側にあった 2 本も落ちた:

| 仮説 | 帰趨 |
|---|---|
| `hfix`(`σ^[p] π = π`) | ★★**落ちた** —— `WildStep.hfix`(§1 の Artin から**無料**) |
| `hiso`(`σ` が等長) | ★★**落ちた**(スペクトルノルムなら) —— `norm_algHom_eq` |

★★**`WildStep` の入力は 9 本、p 進局所体では `hnormp` / `hiso` が自動なので実質 7 本**
(`hi he hdeg htop hjump heM hval`)であり、そこから 10 本すべてを供給して
`_of_conj_pure` を `k = 0` で呼べる(`WildStep.prod_axDecay_via_conj_pure`)。

## ★★④落ちなかったもの(★理由を正確に)

1. ★`hi : 0 < i`(頂点が暴分岐)。★数学的には**全分岐 `p` 次 + 剰余標数 `p` なら自動**だが、
   その導出は「`G_0/G_1 ↪ k^×` かつ `p ∤ |k^×|`」を要し、**剰余体を作る必要がある**
   (本ファイルはノルム言語だけで書かれており剰余体を持たない)。
   ★`RamificationJumpBound` / `GainedJumpSeq` の道にも `i` の**下界**は無い
   (`sub_one_mul_jump_le_of_splits` は**上界** `(p−1)i ≤ p^{k+1}e` のみ)。
2. ★`heM` / `hval` / `hdeg` / `htop` / `he`。★これらは
   `Skeleton/PGC/Setup.lean:40` の `PAdicLocalField`
   (`carrier / Field / Algebra ℚ_p / FiniteDimensional` の**4 つだけ**)からは出ない。
   ★出すには「`M` の中に全分岐巡回 `p` 次の層 `F ⊂ M` と素元 `π` を取る」
   =局所体の構造定理(Eisenstein 多項式・一意化元)が要り、それは本ファイルの持ち場ではない。
   ★★つまり `WildStep` を `PAdicLocalField` から**作る**ことは本波では**できていない**。
   できたのは `hnormp` と `hiso` の 2 本の供給だけである。

## ★★⑤在庫の測定(コマンドつき)

| 部品 | 在庫 | 測定 |
|---|---|---|
| `AlgHom.bijective` | ★mathlib に在った | `grep -n "AlgHom.bijective" .cache/mathlib-index.txt` |
| Artin `[M : M^H] = |H|` | ★mathlib に在った | `IntermediateField.finrank_fixedField_eq_card`(同上) |
| `Nat.card (zpowers a) = orderOf a` | ★mathlib に在った | `Nat.card_zpowers`(同上) |
| `(s.map (X − C ·)).prod.roots = s` | ★mathlib に在った | `Polynomial.roots_multiset_prod_X_sub_C`(同上) |
| スペクトルノルムの Galois 不変性 | ★★**mathlib に在った** | `spectralNorm_eq_of_equiv`(同上) |
| `‖x‖ = spectralNorm k M x` | ★mathlib に在った | `NormedAlgebra.norm_eq_spectralNorm`(同上) |
| 軌道 = 根、跳びの一様性 | ★★木に在った | `grep -n "map_minpoly_eq_prod_iterate\|norm_iterate_sub_self_eq_of_lt_of_fix" lean/ABC3/Found/PGC/MinpolyOrbitSplit.lean` |
| `k = 0` の出口 | ★★木に在った | 同ファイル `exists_norm_sub_algebraMap_le_axDecay_of_orbit` |

★★**「形が嘘」を 1 件踏んだ**(`lean-idioms.md` #297 と同型):

```
error: Function expected at
  IntermediateField.finrank_top'
but this term has type
  finrank ?m.129 ↥⊤ = finrank ?m.129 ?m.131
```

索引の行 `theorem finrank_top' : finrank F (⊤ : IntermediateField F E) = finrank F E` は
section の `variable` を含まないため**明示引数があるように見える**が、実際は引数ゼロである。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. ★`GainedJumpFree` / `GainedJumpSeq` / `MinpolyOrbitSplit` / `LocalFieldNorm` /
   `Skeleton/PGC/Setup.lean` は**読むだけ**で 1 行も書き換えていない。
2. ★構造の場の名前は `π` / `σ` ではなく `unif` / `gen` にした
   (`S.π` は `Real.pi` の記法と紛れうるため)。原典の記号との対応は各 docstring に書いた。
3. ★`WildStep.heM` は `‖p‖ = ‖π‖^{p·e}` と `k = 0` に固定した。
   ★`k ≥ 1` は上の①により空虚なので、一般の `k` を持たせても使い道が無いからである。
-/
namespace ABC3.Found.PGC

namespace PureStepSetup

open Module

/-! ## §1 抽象核(体論だけ) —— `[M : F] = n` なら `τ ^ n = 1` -/

section FieldCore

variable {F M : Type*} [Field F] [Field M] [Algebra F M]

/-- `natDegree (minpoly F π) = n`(`0 < n`)なら `π` は整。 -/
theorem isIntegral_of_natDegree_pos {n : ℕ} (hn : 0 < n) {π : M}
    (hdeg : (minpoly F π).natDegree = n) : IsIntegral F π := by
  by_contra hc
  rw [minpoly.eq_zero hc, Polynomial.natDegree_zero] at hdeg
  omega

theorem finrank_eq_of_adjoin_eq_top {n : ℕ} (hn : 0 < n) {π : M}
    (hdeg : (minpoly F π).natDegree = n)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) : finrank F M = n := by
  have hint : IsIntegral F π := isIntegral_of_natDegree_pos hn hdeg
  have hsub : (IntermediateField.adjoin F ({π} : Set M)).toSubalgebra = ⊤ := by
    rw [IntermediateField.adjoin_simple_toSubalgebra_of_isAlgebraic hint.isAlgebraic]; exact htop
  have h1 : IntermediateField.adjoin F ({π} : Set M) = ⊤ := by
    refine eq_top_iff.mpr fun x _ => ?_
    have hx : x ∈ (IntermediateField.adjoin F ({π} : Set M)).toSubalgebra := by
      rw [hsub]; trivial
    exact hx
  calc finrank F M = finrank F (⊤ : IntermediateField F M) := IntermediateField.finrank_top'.symm
    _ = finrank F (IntermediateField.adjoin F ({π} : Set M)) := by rw [h1]
    _ = n := by rw [IntermediateField.adjoin.finrank hint, hdeg]

theorem orderOf_dvd_finrank (τ : M ≃ₐ[F] M) [FiniteDimensional F M] :
    orderOf τ ∣ finrank F M := by
  refine ⟨finrank F (IntermediateField.fixedField (Subgroup.zpowers τ)), ?_⟩
  rw [← Module.finrank_mul_finrank F (IntermediateField.fixedField (Subgroup.zpowers τ)) M,
    IntermediateField.finrank_fixedField_eq_card, Nat.card_zpowers, mul_comm]

theorem algHom_pow_apply_eq_self_of_natDegree {n : ℕ} (hn : 0 < n) {π : M}
    (hdeg : (minpoly F π).natDegree = n)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) (τ : M →ₐ[F] M) (x : M) :
    (τ ^ n) x = x := by
  have hrank : finrank F M = n := finrank_eq_of_adjoin_eq_top hn hdeg htop
  have hfd : FiniteDimensional F M := FiniteDimensional.of_finrank_pos (by rw [hrank]; exact hn)
  set e : M ≃ₐ[F] M := AlgEquiv.ofBijective τ (AlgHom.bijective τ) with he
  have hdvd : orderOf e ∣ n := by rw [← hrank]; exact orderOf_dvd_finrank e
  have h1 : e ^ n = 1 := orderOf_dvd_iff_pow_eq_one.mp hdvd
  have hcoe : ⇑e = ⇑τ := AlgEquiv.coe_ofBijective τ (AlgHom.bijective τ)
  have h2 : (⇑τ)^[n] x = x := by
    have := congrArg (fun f : M ≃ₐ[F] M => f x) h1
    simpa [AlgEquiv.coe_pow, hcoe] using this
  rw [AlgHom.coe_pow] at *
  exact h2

/-- ★★★**抽象核の系**: `M = F(π)` が `n` 次なら `τ ^ n = 1`。 -/
theorem algHom_pow_eq_one_of_natDegree {n : ℕ} (hn : 0 < n) {π : M}
    (hdeg : (minpoly F π).natDegree = n)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) (τ : M →ₐ[F] M) : τ ^ n = 1 :=
  AlgHom.ext fun y => by
    rw [algHom_pow_apply_eq_self_of_natDegree hn hdeg htop τ y]
    rfl

/-- ★★★★★★**`GainedJumpFree` の `hne` は `1 ≤ k` では偽である**。

`(minpoly F π).natDegree = p` と `Algebra.adjoin F {π} = ⊤` は `[M : F] = p` を意味し、
`τ : M →ₐ[F] M` は(有限次なので)自己同型、その位数は Artin により `p` を割る。
ゆえに `τ ^ p = 1`、したがって `1 ≤ k` では `τ ^ (p ^ k) = 1` で `π` を動かせない。 -/
theorem algHom_pow_pow_apply_eq_self {n k : ℕ} (hn : 0 < n) (hk : 1 ≤ k) {π : M}
    (hdeg : (minpoly F π).natDegree = n)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) (τ : M →ₐ[F] M) (x : M) :
    (τ ^ n ^ k) x = x := by
  have h1 : τ ^ n = 1 := algHom_pow_eq_one_of_natDegree hn hdeg htop τ
  have hpk : n ^ k = n * n ^ (k - 1) := by
    conv_lhs => rw [show k = 1 + (k - 1) by omega]
    rw [pow_add, pow_one]
  rw [hpk, pow_mul, h1, one_pow]
  rfl

end FieldCore

/-! ## §2 `GainedJumpFree` の `_pure` は `1 ≤ k` で空虚 -/

section Vacuous

variable {F M : Type*} [Field F] [NormedField M] [IsUltrametricDist M] [Algebra F M]

omit [IsUltrametricDist M] in
/-- ★★★★★★**`GainedJumpFree.exists_norm_sub_algebraMap_le_prod_axDecay_of_conj_pure` の
仮説は `1 ≤ k` で両立しない。** -/
theorem conj_pure_absurd_of_one_le {p : ℕ} [Fact p.Prime] {k : ℕ} (hk : 1 ≤ k)
    {π : M} (τ : M →ₐ[F] M)
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤)
    (hne : (τ ^ p ^ k) π ≠ π) : False :=
  hne (algHom_pow_pow_apply_eq_self (Fact.out : p.Prime).pos hk hdeg htop τ π)

end Vacuous


/-! ## §3 `hiso`(等長性)の供給 —— スペクトルノルムなら**無料** -/

section Isometry

/-- ★★`hiso` の抽象核: `M` のノルムが完備な基礎体 `k` 上のスペクトルノルムと一致するなら、
`k`-代数同型はすべて等長。★mathlib の `spectralNorm_eq_of_equiv` そのもの。 -/
theorem norm_algEquiv_eq {k M : Type*} [NontriviallyNormedField k] [IsUltrametricDist k]
    [CompleteSpace k] [NormedField M] [NormedAlgebra k M] [Algebra.IsAlgebraic k M]
    (σ : M ≃ₐ[k] M) (y : M) : ‖σ y‖ = ‖y‖ := by
  rw [NormedAlgebra.norm_eq_spectralNorm k (σ y), NormedAlgebra.norm_eq_spectralNorm k y,
    ← spectralNorm_eq_of_equiv σ y]

/-- ★★★`WildStep.hiso` の供給。`σ : M →ₐ[F] M` は `k`-線形でもあり、
`[M : F] < ∞` から全単射なので `M ≃ₐ[k] M` に持ち上がる。 -/
theorem norm_algHom_eq {k F M : Type*} [NontriviallyNormedField k] [IsUltrametricDist k]
    [CompleteSpace k] [Field F] [NormedField M] [Algebra k F] [Algebra F M] [NormedAlgebra k M]
    [IsScalarTower k F M] [Algebra.IsAlgebraic k M] [FiniteDimensional F M]
    (σ : M →ₐ[F] M) (y : M) : ‖σ y‖ = ‖y‖ := by
  have hb : Function.Bijective (σ.restrictScalars k) := AlgHom.bijective σ
  have h := norm_algEquiv_eq (AlgEquiv.ofBijective (σ.restrictScalars k) hb) y
  rwa [AlgEquiv.coe_ofBijective] at h

end Isometry

/-! ## §4 設定を 1 つの構造にまとめる —— `WildStep` -/

/-- ★★★**全分岐巡回 `p` 次の 1 段**の設定。

`GainedJumpFree` の配管 10 本
(`he hnormp heM hval hdeg htop hsep hsp hconj hne`)を 1 つの構造にまとめたもの。
★入力側は 9 本で、うち `hsep` / `hsp` / `hconj` / `hne` は**含まれていない**
——それらは下で**定理として出す**。 -/
structure WildStep (p : ℕ) (F M : Type*) [Field F] [NormedField M] [IsUltrametricDist M]
    [Algebra F M] where
  /-- 素元 -/
  unif : M
  /-- `Gal(M/F)` の生成元 -/
  gen : M →ₐ[F] M
  /-- 跳び -/
  i : ℕ
  /-- `F` の絶対分岐指数 -/
  e : ℕ
  hi : 0 < i
  he : 0 < e
  hdeg : (minpoly F unif).natDegree = p
  htop : Algebra.adjoin F ({unif} : Set M) = ⊤
  hiso : ∀ y : M, ‖gen y‖ = ‖y‖
  hjump : ‖gen unif - unif‖ = ‖unif‖ ^ (i + 1)
  hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹
  heM : ‖(p : M)‖ = ‖unif‖ ^ (p * e)
  hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖unif‖ ^ ((p : ℤ) * m)

namespace WildStep

variable {p : ℕ} [Fact p.Prime] {F M : Type*} [Field F] [NormedField M] [IsUltrametricDist M]
  [Algebra F M] (S : WildStep p F M)

theorem hπ0 : 0 < ‖S.unif‖ := by
  have hp' : p.Prime := Fact.out
  have hpR : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp'.one_lt
  have hEne : p * S.e ≠ 0 := Nat.mul_ne_zero hp'.ne_zero S.he.ne'
  have hπpow : ‖S.unif‖ ^ (p * S.e) = ((p : ℝ))⁻¹ := by rw [← S.heM]; exact S.hnormp
  rcases lt_or_eq_of_le (norm_nonneg S.unif) with h | h
  · exact h
  · exfalso
    rw [← h, zero_pow hEne] at hπpow
    exact absurd hπpow.symm (ne_of_gt (inv_pos.mpr (lt_trans zero_lt_one hpR)))

theorem hπ1 : ‖S.unif‖ < 1 := by
  have hp' : p.Prime := Fact.out
  have hpR : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp'.one_lt
  have hπpow : ‖S.unif‖ ^ (p * S.e) = ((p : ℝ))⁻¹ := by rw [← S.heM]; exact S.hnormp
  refine GainedBridgeSupply.lt_one_of_pow_lt_one (E := p * S.e) ?_
  rw [hπpow, inv_lt_one_iff₀]; exact Or.inr hpR

theorem hint : IsIntegral F S.unif :=
  isIntegral_of_natDegree_pos (Fact.out : p.Prime).pos S.hdeg

/-- ★★`σ^[p] π = π`(`GainedJumpFree` では**仮説**だったが、
`[M : F] = p` と Artin から**無料で出る**)。 -/
theorem hfix : (S.gen : M → M)^[p] S.unif = S.unif := by
  have h := algHom_pow_apply_eq_self_of_natDegree (Fact.out : p.Prime).pos S.hdeg S.htop
    S.gen S.unif
  rwa [AlgHom.coe_pow] at h

theorem hne : S.gen S.unif ≠ S.unif := by
  intro h
  have hj := S.hjump
  rw [h, sub_self, norm_zero] at hj
  exact absurd hj.symm (ne_of_gt (pow_pos S.hπ0 _))

/-- ★★`GainedJumpFree` の `hsp`。★構造の入力には**入っていない**。 -/
theorem hsp : ((minpoly F S.unif).map (algebraMap F M)).Splits := by
  have hp' : p.Prime := Fact.out
  refine Polynomial.splits_iff_card_roots.mpr ?_
  rw [(minpoly.monic S.hint).natDegree_map, S.hdeg,
    map_minpoly_eq_prod_iterate hp' S.gen S.hfix S.hne S.hdeg,
    Polynomial.roots_multiset_prod_X_sub_C]
  simp

/-- ★★`GainedJumpFree` の `hsep`。★構造の入力には**入っていない**。 -/
theorem hsep : (minpoly F S.unif).Separable := by
  have hp' : p.Prime := Fact.out
  rw [← Polynomial.separable_map (algebraMap F M),
    map_minpoly_eq_prod_iterate hp' S.gen S.hfix S.hne S.hdeg]
  have hprod : (Multiset.map (fun a => Polynomial.X - Polynomial.C a)
      ((Multiset.range p).map fun k => (S.gen : M → M)^[k] S.unif)).prod
      = ∏ k ∈ Finset.range p,
          (Polynomial.X - Polynomial.C ((S.gen : M → M)^[k] S.unif)) := by
    rw [Multiset.map_map]; rfl
  rw [hprod]
  refine Polynomial.separable_prod_X_sub_C_iff'.mpr fun x hx y hy hxy => ?_
  exact injOn_iterate_of_prime hp' S.hfix S.hne (Finset.mem_range.mp hx)
    (Finset.mem_range.mp hy) hxy

/-- ★★★`GainedJumpFree` の `hconj`(全共役が同じ距離)。★構造の入力には**入っていない**。 -/
theorem hconj (a : M) (ha : ((minpoly F S.unif).map (algebraMap F M)).IsRoot a)
    (hane : a ≠ S.unif) : ‖S.unif - a‖ = ‖S.unif‖ ^ (S.i + 1) := by
  have hp' : p.Prime := Fact.out
  have hmonic : ((minpoly F S.unif).map (algebraMap F M)).Monic :=
    (minpoly.monic S.hint).map _
  have hmem : a ∈ ((minpoly F S.unif).map (algebraMap F M)).roots :=
    Polynomial.mem_roots'.mpr ⟨hmonic.ne_zero, ha⟩
  rw [map_minpoly_eq_prod_iterate hp' S.gen S.hfix S.hne S.hdeg,
    Polynomial.roots_multiset_prod_X_sub_C] at hmem
  obtain ⟨k, hk, rfl⟩ := Multiset.mem_map.mp hmem
  have hkp : k < p := Multiset.mem_range.mp hk
  have hk0 : 0 < k := by
    rcases Nat.eq_zero_or_pos k with rfl | h
    · exact absurd (by simp) hane
    · exact h
  rw [norm_sub_rev]
  exact norm_iterate_sub_self_eq_of_lt_of_fix (fun x y => map_sub S.gen x y) S.hiso hp'
    S.hfix hk0 hkp S.hjump

/-! ### 出口 -/

/-- ★★★★**構造 1 つから `axDecay p 1` の 1 段が出る**。 -/
theorem exists_norm_sub_algebraMap_le_axDecay (x : M) :
    ∃ y : F, ‖x - algebraMap F M y‖ ≤ axDecay p 1 * ‖S.gen x - x‖ := by
  rw [axDecay_one]
  exact exists_norm_sub_algebraMap_le_axDecay_of_orbit (Fact.out : p.Prime) S.gen S.hπ0 S.hπ1
    S.hi S.hiso S.hfix S.hjump S.hdeg S.htop S.hval S.hnormp x

/-- ★★★★★**`GainedJumpFree...of_conj_pure` の結論の形**(★`k = 0`)。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay (x : M) :
    ∃ y : F, ‖x - algebraMap F M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (0 + 1), axDecay p j) * ‖S.gen x - x‖ := by
  simpa using S.exists_norm_sub_algebraMap_le_axDecay x

/-- ★★★★★★**配られた配管 10 本を構造 1 つから供給して
`GainedJumpFree.exists_norm_sub_algebraMap_le_prod_axDecay_of_conj_pure` を呼ぶ**(`k = 0`)。

★`he hnormp heM hval hdeg htop` は構造の場、
★`hsep hsp hconj hne` は上で証明した定理である。
★★`k ≥ 1` では `hne` が偽なので(`conj_pure_absurd_of_one_le`)、
★★これが `_of_conj_pure` を非空虚に使える**唯一の**場合である。 -/
theorem prod_axDecay_via_conj_pure (x : M) :
    ∃ y : F, ‖x - algebraMap F M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (0 + 1), axDecay p j) * ‖S.gen x - x‖ :=
  GainedJumpFree.exists_norm_sub_algebraMap_le_prod_axDecay_of_conj_pure
    (e := S.e) (i := S.i) (k := 0) S.gen S.he S.hi S.hnormp (by simpa using S.heM)
    S.hval S.hdeg S.htop S.hsep S.hsp S.hconj (by simpa using S.hne) x

end WildStep

/-! ## §5 `PAdicLocalField` からの供給 -/

section PAdic

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-- ★★`WildStep.hnormp` の供給。★`PAdicLocalField` が持つ 4 つ
(`carrier / Field / Algebra ℚ_p / FiniteDimensional`)だけから出る。

★`LocalFieldNorm.norm_natCast_p_lt_one` は `< 1` までしか言っていなかったが、
その証明はそのまま**等式**を出す。 -/
theorem norm_natCast_p_eq_inv (K : PAdicLocalField p) :
    ‖((p : ℕ) : K.carrier)‖ = ((p : ℝ))⁻¹ := by
  rw [show ((p : ℕ) : K.carrier) = algebraMap ℚ_[p] K.carrier ((p : ℕ) : ℚ_[p]) from
        (map_natCast _ p).symm,
      norm_algebraMap, Padic.norm_p]

/-- ★★`WildStep.hiso` の供給(`PAdicLocalField` 版)。
`F` が `ℚ_p` と `K.carrier` の間の中間層なら、`K.carrier →ₐ[F] K.carrier` は等長。 -/
theorem norm_algHom_eq_of_pAdicLocalField (K : PAdicLocalField p) {F : Type*} [Field F]
    [Algebra ℚ_[p] F] [Algebra F K.carrier] [IsScalarTower ℚ_[p] F K.carrier]
    [FiniteDimensional F K.carrier] (σ : K.carrier →ₐ[F] K.carrier) (y : K.carrier) :
    ‖σ y‖ = ‖y‖ :=
  norm_algHom_eq (k := ℚ_[p]) σ y

end PAdic

/-! ## §6 `.src`(原典の対応箇所) -/

def algHom_pow_pow_apply_eq_self.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def conj_pure_absurd_of_one_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def WildStep.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def WildStep.prod_axDecay_via_conj_pure.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §7 使っている公理の一覧 -/

#print axioms isIntegral_of_natDegree_pos
#print axioms finrank_eq_of_adjoin_eq_top
#print axioms orderOf_dvd_finrank
#print axioms algHom_pow_apply_eq_self_of_natDegree
#print axioms algHom_pow_eq_one_of_natDegree
#print axioms algHom_pow_pow_apply_eq_self
#print axioms conj_pure_absurd_of_one_le
#print axioms norm_algEquiv_eq
#print axioms norm_algHom_eq
#print axioms WildStep.hfix
#print axioms WildStep.hne
#print axioms WildStep.hsp
#print axioms WildStep.hsep
#print axioms WildStep.hconj
#print axioms WildStep.exists_norm_sub_algebraMap_le_axDecay
#print axioms WildStep.prod_axDecay_via_conj_pure
#print axioms norm_natCast_p_eq_inv
#print axioms norm_algHom_eq_of_pAdicLocalField

end PureStepSetup

end ABC3.Found.PGC
