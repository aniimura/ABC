import ABC3.Found.PGC.LubinTateUpperRamificationVanish
import ABC3.Found.PGC.LubinTateReciprocityIsomorphism
import ABC3.Found.PGC.LubinTateFormalGroupLawEstimate
import ABC3.Found.PGC.LubinTateTowerCompatible
import ABC3.Found.PGC.LubinTatePsiNorm
import ABC3.Found.PGC.LubinTateDegree
import ABC3.Found.PGC.UnramifiedBaseChangeInvariance

/-!
# `i(σ) = q^{v_K(u−1)}` —— Yoshida 2008 Proposition 6.14 の証明の段 1

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Proposition 6.14(物理 p.17)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
の `section-6.html` の `#prop-6-14`。

原文 (Yoshida08 p.17):
> Proposition 6.14. Let x ∈ K^× with v(x) = n > 0. Let L = K_n and K^m_x as in Definition 5.3. Then we have Gal(K^m_x/L)^m = {id} for all m ≥ 1 (see Definition 6.12).

`Found/PGC/LubinTateUpperRamificationVanish.lean` が Prop 6.14(`n = 1` 形)を
**段 1 だけを仮定 `hρ` に残して**閉じている。本ファイルはその段 1 を埋める。

## 原典の段 1(全 4 行。★`>` を付けない ―― 逐語引用は上の 1 本だけである)

    For σ ≠ id, set β := [u−1]_f(α). If v_K(u−1) = i for 0 ≤ i < m, then β ∈ µ^×_{f,m−i}
    by Lemma 4.3(ii). Hence β is a uniformizer of K^{m−i}_x by Proposition 4.4(ii), which
    shows v(β) = q^i. Now σ(α) = [u]_f(α) = α +_f β ≡ α+β (mod αβ), hence
    i(σ) = v(σ(α)−α) = v(β) = q^i.

## 何が入ったか

原典の段 1 は 3 つに割れる。**3 つとも入った。**

| # | 主張 | 本ファイル |
|---|---|---|
| 1a | `v_K(a) = i` ⟹ `[a]_f(α) ∈ µ^×_{f,m−i}` | `action_pi_pow_mul_unit_mem_iteratedLubinTatePsiTorsionPoints` |
| 1b | `v(β) = q^i` | `norm_smul_torsionGen_sub`(★**塔の分岐指数を経由しない**、下記) |
| 1c | `F_f(X,Y) ≡ X+Y (mod XY)` から `v(σα−α) = v(β)` | `norm_sub_eq_of_norm_sub_sub_le_mul` + `norm_lubinTateAction_one_add_sub` |

到達点は `ramIndex_torsionGen_eq_pow`:

    `i(σ) = q^i`   (`σ(α) = [u]_f(α)`, `u − 1 = π^i·w`, `w ∈ 𝒪_K^×`, `m = i+k+1`)

★**等号である**(`≥` ではない)。★**下付き** `i(σ)`(`ramIndex`)の話であって
上付きではない。

## 1c は「壁」ではなかった —— 木に既にあった

段 1 を割った実装者は 1c を「壁」と見立て、「木の `formalGroupLaw` は 1 次係数しか
押さえておらず、2 変数冪級数の『混合項は `XY` で割れる』補題が無い」と報告した。
★**測ったら在った** ―― `Found/PGC/LubinTateFormalGroupLawEstimate.lean` の

    `norm_aeval_formalGroupLaw_sub_le : ‖F_f(z,w) − z − w‖ ≤ ‖z‖·‖w‖`

が、まさに「混合項は `XY` で割れる」の**評価版**である(`F_f(X,0)=X` と `F_f(0,Y)=Y`
の**両方**を使って `i≥1 ∧ j≥1` を出しており、片方だけでは出ない)。
★**木の `formalGroupLaw` は `MvPowerSeries (Fin 2) A`**(`X 0`・`X 1` の 2 変数)で
あり、係数は `Finsupp.single (0 : Fin 2) n` などで取る。

残っていたのは**評価不等式から等号 `v(σα−α) = v(β)` を出す 15 行**だけで、
それが本ファイルの抽象核 `norm_sub_eq_of_norm_sub_sub_le_mul` である。

## 抽象核(§0。分岐・付値・Lubin-Tate・Galois の語彙が 1 つも出てこない)

* `norm_sub_eq_of_norm_sub_sub_le_mul`
  —— 超距離ノルム環で `‖s − a − b‖ ≤ ‖a‖·‖b‖` かつ `‖a‖ < 1` ならば `‖s − a‖ = ‖b‖`。
  ★`b = 0` でも成り立つので `b ≠ 0` は要らない(★退化していない: `‖a‖ < 1` が本質で、
  これを外すと `‖a‖ = 1` のとき `‖s−a‖ = ‖b‖` は偽になり得る)。
* `addVal_eq_of_norm_eq_pow`
  —— DVR `B` 上の乗法的で `≤ 1` の実数値 `N` について、`N y = N ϖ ^ n`(`ϖ` 既約、
  `0 < N ϖ < 1`)ならば `addVal B y = n`。★**ノルム語と付値語の橋**。
* `not_pow_succ_dvd_pow_mul_unit` —— 整域で `ϖ^{k+1} ∤ ϖ^k·(単数)`。
* `pow_rpow_one_div_mul` —— `(c^{1/(N·d)})^N = c^{1/d}`(実数の `rpow` だけ)。
* `nat_pow_sub_factor` —— `q^{i+k+1} − q^{i+k} = q^i·(q^{k+1} − q^k)`(ℕ の等式)。

## 原典より短い道 —— 段 1b は塔の分岐指数を経由しない

原典は「`β` は `K^{m−i}_x` の素元だから `v(β) = q^i`」と、**塔 `K^m_x/K^{m−i}_x` の
分岐指数が `q^i`** であることを使う。★本ファイルはそれを使わない ―― `α` も `β` も
**`K` 上のスペクトルノルム**で測り、

    `‖α‖ = ‖π‖^{1/deg ψ_m}`,  `‖β‖ = ‖π‖^{1/deg ψ_{m−i}}`,  `deg ψ_n = q^n − q^{n−1}`

から `‖β‖ = ‖α‖^{q^i}` を出すだけで済む(`spectralNorm_root_iteratedLubinTatePsi` +
`natDegree_iteratedLubinTatePsi`)。★**塔も相対 Lubin-Tate も要らない。**

## 逸脱の記録

1. ★**`ρ` による翻訳 `σ(α) = [u]_f(α)` は仮定ではなく `galoisUnitReciprocityMap_spec`
   から出している**(`smul_torsionGen_eq_lubinTateAction`)。原典の Proposition 4.4(iii)
   に対応する。
2. ★**`u − 1 = π^i · w`(`w ∈ 𝒪_K^×`)という形で `v_K(u−1) = i` を表している。**
   原典の `v_K(u−1) = i` と同値だが、ℕ の切り詰め引き算・除算を出さないための書き換え。
   `0 ≤ i < m` は `m = i + k + 1` と書くことで埋め込んである。
3. ★**`α` は `x` 自身**(`torsionGen`)である。原典の「`α ∈ µ^×_{f,m}` は `K^m_x` の
   素元」は木に在庫が無かった(`isTotallyRamifiedAdjoin_iteratedLubinTatePsi` は
   完全分岐までしか言わない)ので、本ファイルが `irreducible_torsionGen` として
   示した。★**仮定ではない。**
4. ★`m = 0` は除いてある(`m = i+k+1 ≥ 1`)。`m = 0` では `µ_{f,0} = {0}` で拡大が
   自明になり、主張が空虚になる。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-! ### §0 抽象核

分岐・付値・Lubin-Tate・Galois の語彙が 1 つも出てこない 5 本。 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**抽象核(1c の心臓)**——超距離ノルム環で
`‖s − a − b‖ ≤ ‖a‖·‖b‖` かつ `‖a‖ < 1` ならば `‖s − a‖ = ‖b‖`。

原典の「`σ(α) = α +_f β ≡ α+β (mod αβ)` ゆえ `v(σα−α) = v(β)`」がこれである。
★`‖a‖ < 1` を落とすと偽になり得る(`‖a‖ = 1` なら `‖a‖·‖b‖ = ‖b‖` で、
超距離の等号が壊れる)。★`b = 0` の場合も込みで成り立つ。 -/
theorem norm_sub_eq_of_norm_sub_sub_le_mul {S : Type*} [NormedCommRing S] [IsUltrametricDist S]
    (s a b : S) (h : ‖s - a - b‖ ≤ ‖a‖ * ‖b‖) (ha : ‖a‖ < 1) :
    ‖s - a‖ = ‖b‖ := by
  rcases eq_or_ne b 0 with rfl | hb
  · simp only [norm_zero, mul_zero, sub_zero, norm_le_zero_iff] at h
    simp [h]
  · have hbpos : 0 < ‖b‖ := norm_pos_iff.mpr hb
    have hlt : ‖s - a - b‖ < ‖b‖ := by
      calc ‖s - a - b‖ ≤ ‖a‖ * ‖b‖ := h
        _ < 1 * ‖b‖ := by exact mul_lt_mul_of_pos_right ha hbpos
        _ = ‖b‖ := one_mul _
    have hne : ‖b‖ ≠ ‖s - a - b‖ := ne_of_gt hlt
    have hmax := IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (x := b) (y := s - a - b) hne
    have hsum : b + (s - a - b) = s - a := by ring
    rw [hsum] at hmax
    rw [hmax, max_eq_left hlt.le]

/-- **抽象核**——整域で `ϖ` が非単数かつ `ϖ ≠ 0` なら `ϖ^{k+1} ∤ ϖ^k·w`(`w` は単数)。 -/
theorem not_pow_succ_dvd_pow_mul_unit {A : Type*} [CommRing A] [IsDomain A] {ϖ : A} (hϖ0 : ϖ ≠ 0)
    (hϖu : ¬ IsUnit ϖ) (k : ℕ) (w : Aˣ) : ¬ (ϖ ^ (k + 1) ∣ ϖ ^ k * (w : A)) := by
  intro hdvd
  obtain ⟨c, hc⟩ := hdvd
  have hcancel : ϖ ^ k * (w : A) = ϖ ^ k * (ϖ * c) := by rw [hc]; ring
  have h1 : (w : A) = ϖ * c := mul_left_cancel₀ (pow_ne_zero k hϖ0) hcancel
  have h2 : ϖ * (c * (↑w⁻¹ : A)) = 1 := by
    rw [← mul_assoc, ← h1, Units.mul_inv]
  exact hϖu (isUnit_iff_exists_inv.mpr ⟨_, h2⟩)

/-- **抽象核の補助**——乗法的で `≤ 1` の `N` について `N(u·ϖ^j) = (N ϖ)^j`(`u` は単数)。 -/
theorem norm_unit_mul_pow_eq {B : Type*} [CommRing B]
    {N : B → ℝ} (hnonneg : ∀ a : B, 0 ≤ N a)
    (hmul : ∀ a b : B, N (a * b) = N a * N b) (hle : ∀ a : B, N a ≤ 1)
    {ϖ : B} (hϖpos : 0 < N ϖ) (u : Bˣ) (j : ℕ) : N ((u : B) * ϖ ^ j) = N ϖ ^ j := by
  have hN1 : N 1 = 1 := by
    have h := hmul ϖ 1
    rw [mul_one] at h
    nlinarith [h]
  have hNu : N (u : B) = 1 := by
    have h1 : N (u : B) * N ((u⁻¹ : Bˣ) : B) = 1 := by
      rw [← hmul]; simpa using hN1
    nlinarith [hnonneg (u : B), hnonneg ((u⁻¹ : Bˣ) : B), hle (u : B), hle ((u⁻¹ : Bˣ) : B)]
  rw [hmul, hNu, one_mul]
  induction j with
  | zero => simpa using hN1
  | succ t ih => rw [pow_succ, hmul, ih, pow_succ]

/-- ★★★★★★★★★★**抽象核——付値語からノルム語へ**: `addVal B y = n` ならば
`N y = (N ϖ)^n`。★こちら向きは `N ϖ < 1` を要らない。 -/
theorem norm_eq_pow_of_addVal_eq {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B]
    {N : B → ℝ} (hnonneg : ∀ a : B, 0 ≤ N a)
    (hmul : ∀ a b : B, N (a * b) = N a * N b) (hle : ∀ a : B, N a ≤ 1)
    {ϖ : B} (hϖ : Irreducible ϖ) (hϖpos : 0 < N ϖ)
    {y : B} {n : ℕ} (h : IsDiscreteValuationRing.addVal B y = (n : ℕ∞)) :
    N y = N ϖ ^ n := by
  rcases eq_or_ne y 0 with rfl | hy0
  · rw [IsDiscreteValuationRing.addVal_zero] at h
    exact absurd h.symm (ENat.coe_ne_top n)
  · obtain ⟨j, u, hju⟩ := IsDiscreteValuationRing.eq_unit_mul_pow_irreducible hy0 hϖ
    rw [IsDiscreteValuationRing.addVal_def y u hϖ j hju] at h
    have hjn : j = n := by exact_mod_cast h
    rw [hju, norm_unit_mul_pow_eq hnonneg hmul hle hϖpos, hjn]

/-- ★★★★★★★★★★★★★★★★★★★★**抽象核——ノルム語から付値語へ**: DVR `B` 上の実数値関数
`N` が `N 0 = 0`・非負・乗法的・`≤ 1` を満たし `ϖ` が既約で `0 < N ϖ < 1` のとき、
`N y = (N ϖ)^n` ならば `addVal B y = n`。

★木には「スペクトルノルム ↔ `addVal`」の橋が無かったので新しく入れた。
★`N` を抽象化してあるので、`adjoinIntegers K x` のノルムに限らず使える。
★★**本ファイルは両向きとも使う** —— こちらは `i(σ) = q^i` の結論を出すのに、
`norm_eq_pow_of_addVal_eq` は `‖π‖ = ‖ϖ‖^e`(素元の分岐指数)を出すのに使う。 -/
theorem addVal_eq_of_norm_eq_pow {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B]
    {N : B → ℝ} (hzero : N 0 = 0) (hnonneg : ∀ a : B, 0 ≤ N a)
    (hmul : ∀ a b : B, N (a * b) = N a * N b) (hle : ∀ a : B, N a ≤ 1)
    {ϖ : B} (hϖ : Irreducible ϖ) (hϖpos : 0 < N ϖ) (hϖlt : N ϖ < 1)
    {y : B} {n : ℕ} (h : N y = N ϖ ^ n) :
    IsDiscreteValuationRing.addVal B y = (n : ℕ∞) := by
  have hy0 : y ≠ 0 := by
    intro hy
    rw [hy, hzero] at h
    exact absurd h.symm (ne_of_gt (pow_pos hϖpos n))
  obtain ⟨j, u, hju⟩ := IsDiscreteValuationRing.eq_unit_mul_pow_irreducible hy0 hϖ
  have hNy : N y = N ϖ ^ j := by rw [hju, norm_unit_mul_pow_eq hnonneg hmul hle hϖpos]
  have hpow : N ϖ ^ j = N ϖ ^ n := by rw [← hNy, h]
  have hjn : j = n := (pow_right_strictAnti₀ hϖpos hϖlt).injective hpow
  rw [IsDiscreteValuationRing.addVal_def y u hϖ j hju, hjn]

/-- **抽象核**——実数の `rpow` だけ: `(c^e)^{1/e} = c`。 -/
theorem rpow_one_div_natCast_pow {c : ℝ} (hc : 0 ≤ c) (e : ℕ) (he : e ≠ 0) :
    (c ^ e : ℝ) ^ (1 / (e : ℝ)) = c := by
  have he' : (e : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr he
  rw [← Real.rpow_natCast c e, ← Real.rpow_mul hc, mul_one_div, div_self he', Real.rpow_one]

/-- **抽象核**——実数の `rpow` だけ: `(c^{1/(N·d)})^N = c^{1/d}`。 -/
theorem pow_rpow_one_div_mul {c : ℝ} (hc : 0 ≤ c) (N d : ℕ) (hN : N ≠ 0) (hd : d ≠ 0) :
    (c ^ (1 / ((N * d : ℕ) : ℝ))) ^ N = c ^ (1 / (d : ℝ)) := by
  rw [← Real.rpow_natCast (c ^ (1 / ((N * d : ℕ) : ℝ))) N, ← Real.rpow_mul hc]
  congr 1
  have hN' : (N : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hN
  have hd' : (d : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hd
  push_cast
  field_simp

/-- **抽象核**——ℕ の等式 `q^{i+k+1} − q^{i+k} = q^i·(q^{k+1} − q^k)`。
★切り詰め引き算は `Nat.mul_sub` が分配してくれるので `q ≥ 2` すら要らない。 -/
theorem nat_pow_sub_factor (q i k : ℕ) :
    q ^ (i + k + 1) - q ^ (i + k) = q ^ i * (q ^ (k + 1) - q ^ k) := by
  rw [Nat.mul_sub, ← pow_add, ← pow_add, ← Nat.add_assoc]

/-! ### §1 具体層 —— 1c(混合項) -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**段 1c**: `‖[1+a]_f(α) − [1]_f(α)‖ = ‖[a]_f(α)‖`。

原典の「`σ(α) = [u]_f(α) = α +_f β ≡ α+β (mod αβ)` ゆえ `v(σα−α) = v(β)`」そのもの。
`lubinTateAction_add`(`(1+a)·α = F_f(1·α, a·α)`)と
`norm_aeval_formalGroupLaw_sub_le`(`‖F_f(z,w) − z − w‖ ≤ ‖z‖·‖w‖`)を合わせ、
抽象核 `norm_sub_eq_of_norm_sub_sub_le_mul` に流し込むだけ。
`‖1·α‖ < 1` は `1·α` が捩れ点であること(`spectralNorm_lt_one_of_mem_…`)から。 -/
theorem norm_lubinTateAction_one_add_sub
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (x : K.closure)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (a : 𝒪[K.carrier]) :
    ‖lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (1 + a) -
        lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem 1‖ =
      ‖lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem a‖ := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  have hadd : lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem (1 + a) =
      MvPowerSeries.aeval (hasEval_actionFam2 K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem 1 a)
        (formalGroupLaw hq hπmax f hf0 hf1 hf) :=
    lubinTateAction_add K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem 1 a
  set α := lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem 1 with hα
  set β := lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem a with hβ
  have hbound : ‖MvPowerSeries.aeval
        (hasEval_actionFam2 K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem 1 a)
        (formalGroupLaw hq hπmax f hf0 hf1 hf) -
        (fun i : Fin 2 => if i = 0 then α else β) 0 -
        (fun i : Fin 2 => if i = 0 then α else β) 1‖ ≤
      ‖(fun i : Fin 2 => if i = 0 then α else β) 0‖ *
        ‖(fun i : Fin 2 => if i = 0 then α else β) 1‖ :=
    norm_aeval_formalGroupLaw_sub_le hq hπmax hπne0 f hf0 hf1 hf
      (hasEval_actionFam2 K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem 1 a)
      α.2 β.2 (algebraMap_mem_adjoinIntegers K x)
  simp only [if_pos, if_neg (by decide : (1 : Fin 2) ≠ 0)] at hbound
  rw [← hadd] at hbound
  have hαlt1 : ‖α‖ < 1 := by
    show spectralNorm K.carrier K.closure
      (↑(↑(α : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)) < 1
    exact spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n _
      (lubinTateActionAtTorsionPoint_mem K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem 1)
  exact norm_sub_eq_of_norm_sub_sub_le_mul _ _ _ hbound hαlt1

/-! ### §2 具体層 —— 1a(原始性が `v_K(a)` だけ下がる) -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**段 1a**(原典の「`β ∈ µ^×_{f,m−i}` by Lemma 4.3(ii)」):
`x` が原始的な `π^{i+k+1}`-捩れ点(`ψ_{i+k+1}` の根)で `w ∈ 𝒪_K^×` のとき、
`[π^i·w]_f(x)` は原始的な `π^{k+1}`-捩れ点(`ψ_{k+1}` の根)。

★木には `[π](α)` の 1 段版(`lubinTateActionAtTorsionPoint_pi_mem_iteratedLubinTatePsi
TorsionPoints`)と単数版(`unit_action_mem_iteratedLubinTatePsiTorsionPoints`)しか
無く、この 2 つを合成しようとすると「`[π](α)` 自身の `adjoinIntegers`」が要って
座標系をまたぐ。★**本補題は `x` 自身の座標系だけで閉じる**——上の 2 本と同じ骨格を
`π^i·w` に対して 1 回だけ回す。

* `∈ Λ_{k+1}`: `[π^{k+1}]([π^i w](x)) = [π^{i+k+1}·w](x) = 0`(`π^{i+k+1} ∣ π^{i+k+1}w`)
  から `D_{k+1}([π^i w](x)) = 0`。
* `∉ Λ_k`: もし `D_k([π^i w](x)) = 0` なら `[π^k]_f = D_k·U_k` により
  `[π^{i+k}·w](x) = 0`、すなわち `π^{i+k+1} ∣ π^{i+k}·w` で
  `not_pow_succ_dvd_pow_mul_unit` に矛盾。 -/
theorem action_pi_pow_mul_unit_mem_iteratedLubinTatePsiTorsionPoints
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (i k : ℕ) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1)
      (by omega))
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1))
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (w : (𝒪[K.carrier])ˣ) :
    (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem
        (π ^ i * (w : 𝒪[K.carrier]))) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure) ∈
      iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (k + 1) (by omega) := by
  haveI := completeSpace_adjoinIntegers K x
  haveI := isLinearTopology_adjoinIntegers K x
  haveI := continuousSMul_adjoinIntegers K x
  have hπnu : ¬ IsUnit π := by
    have hmemπ : π ∈ IsLocalRing.maximalIdeal (𝒪[K.carrier]) := by
      rw [hπmax]; exact Ideal.mem_span_singleton_self π
    exact mem_nonunits_iff.mp hmemπ
  set a : 𝒪[K.carrier] := π ^ i * (w : 𝒪[K.carrier]) with ha_def
  set z := lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem a
    with hz_def
  set g : adjoinIntegers K x →+* K.closure :=
    (algebraMap (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure).comp
      (algebraMap (adjoinIntegers K x) (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    with hg_def
  have hginj : Function.Injective g := fun s t h => Subtype.ext (Subtype.ext h)
  have hgcomp : g.comp (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) =
      algebraMap (𝒪[K.carrier]) K.closure := by
    apply RingHom.ext; intro y; rfl
  have hzero : lubinTateEvalAtPoint K x z
      (hasEval_lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem a)
      (LubinTateAction hq hπmax f hf0 hf1 hf (π ^ (k + 1))) = 0 := by
    rw [lubinTateAction_mul K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem (π ^ (k + 1)) a,
      lubinTateActionAtTorsionPoint_eq_zero_iff_dvd_of_mem_iteratedLubinTatePsiTorsionPoints
        K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) (by omega) x hxψ hmem hxn]
    exact ⟨(w : 𝒪[K.carrier]), by rw [ha_def]; ring⟩
  have hDk1 : Polynomial.aeval z
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (k + 1)) = 0 :=
    eq_zero_of_pi_pow_action_eq_zero K hq hπmax hπne0 f hf0 hf1 hf (k + 1) x z _ hzero
  have hzmem : (g z : K.closure) ∈
      iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (k + 1) := by
    rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots']
    refine ⟨(isDistinguishedAt_iteratedLubinTateDistinguished
      hq hπmax hπne0 f hf0 hf1 hf (k + 1)).monic.map _ |>.ne_zero, ?_⟩
    have key := Polynomial.hom_eval₂
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (k + 1))
      (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) g z
    rw [← Polynomial.aeval_def, hgcomp, hDk1, map_zero] at key
    show (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf (k + 1))).eval (g z) = 0
    rw [Polynomial.eval_map]
    exact key.symm
  rw [← iteratedLubinTateTorsionPoints_sdiff_eq_iteratedLubinTatePsiTorsionPoints
    K hq hπmax hπne0 f hf0 hf1 hf (k + 1) (by omega), Finset.mem_sdiff]
  refine ⟨hzmem, fun hbad => ?_⟩
  rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots'] at hbad
  obtain ⟨_, hbadroot⟩ := hbad
  rw [Polynomial.IsRoot, Polynomial.eval_map] at hbadroot
  have hDk : Polynomial.aeval z
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf k) = 0 := by
    apply hginj
    have key := Polynomial.hom_eval₂
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf k)
      (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x)) g z
    rw [← Polynomial.aeval_def, hgcomp] at key
    rw [key, map_zero]
    exact hbadroot
  have hzero' : lubinTateEvalAtPoint K x z
      (hasEval_lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem a)
      (LubinTateAction hq hπmax f hf0 hf1 hf (π ^ k)) = 0 := by
    show PowerSeries.aeval _ (LubinTateAction hq hπmax f hf0 hf1 hf (π ^ k)) = 0
    rw [LubinTateAction_pi_pow hq hπmax hπne0 f hf0 hf1 hf k,
      iteratedLubinTate_eq_distinguished_mul_unit hq hπmax hπne0 f hf0 hf1 hf k,
      map_mul, PowerSeries.aeval_coe, hDk, zero_mul]
  have hprod0 : lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem
      (π ^ k * a) = 0 :=
    (lubinTateAction_mul K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem (π ^ k) a).symm.trans
      hzero'
  rw [lubinTateActionAtTorsionPoint_eq_zero_iff_dvd_of_mem_iteratedLubinTatePsiTorsionPoints
    K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) (by omega) x hxψ hmem hxn] at hprod0
  refine not_pow_succ_dvd_pow_mul_unit hπne0 hπnu (i + k) w ?_
  have hrw : π ^ k * a = π ^ (i + k) * (w : 𝒪[K.carrier]) := by rw [ha_def]; ring
  rwa [hrw] at hprod0

/-! ### §3 具体層 —— 1b(`ψ_n` の根のノルム。★塔の分岐指数を経由しない) -/

/-- **`ψ_n` の根のスペクトルノルムは `‖π‖^{1/(q^n − q^{n−1})}`**——
`spectralNorm_root_iteratedLubinTatePsi`(任意の根について)と
`natDegree_iteratedLubinTatePsi`(`deg ψ_n = q^n − q^{n−1}`)を合わせただけ。
`iteratedLubinTatePsiTorsionPoints` は `ψ_n` を `K.closure` へ写した多項式の根の
`Finset` なので、`Polynomial.mem_roots'` から根であることを取り出して
`K.carrier` 経由の `aeval` へ移す(`Polynomial.eval₂_map` + `algebraMap` の合成)。 -/
theorem spectralNorm_of_mem_iteratedLubinTatePsiTorsionPoints
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (y : K.closure)
    (hy : y ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn) :
    spectralNorm K.carrier K.closure y =
      ‖π‖ ^ (1 / (((pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) : ℕ) : ℝ)) := by
  rw [iteratedLubinTatePsiTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots'] at hy
  obtain ⟨_, hroot⟩ := hy
  rw [Polynomial.IsRoot, Polynomial.eval_map] at hroot
  have hcomp : (algebraMap K.carrier K.closure).comp (algebraMap (𝒪[K.carrier]) K.carrier) =
      algebraMap (𝒪[K.carrier]) K.closure := by
    apply RingHom.ext; intro c; rfl
  have haeval : Polynomial.aeval y (Polynomial.map (algebraMap (𝒪[K.carrier]) K.carrier)
      (iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn)) = 0 := by
    rw [Polynomial.aeval_def, Polynomial.eval₂_map, hcomp]
    exact hroot
  rw [spectralNorm_root_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf n hn y haeval,
    natDegree_iteratedLubinTatePsi hq hπmax hπne0 f hf0 hf1 hf n hn]

/-! ### §4 具体層 —— 段 1 の組み立て -/

/-- 捩れ点 `x` を `adjoinIntegers K x` の元として見たもの(原典の `α ∈ µ^×_{f,m}`)。 -/
noncomputable def torsionGen
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (x : K.closure)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    adjoinIntegers K x :=
  ⟨⟨x, hmem⟩,
    mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints
      K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem⟩

/-- **`ρ` の翻訳**(原典の Proposition 4.4(iii)): `ρ_{f,m}(σ) = u mod p^m` なら
`σ(α) = [u]_f(α)`。`galoisUnitReciprocityMap_spec` と `unitActionQuotientLift_mk` から。 -/
theorem smul_torsionGen_eq_lubinTateAction
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (σ : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≃ₐ[K.carrier]
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (u : (𝒪[K.carrier])ˣ)
    (hu : galoisUnitReciprocityMap K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ =
      (u : (𝒪[K.carrier])ˣ ⧸ principalUnits K π n)) :
    σ • torsionGen K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem =
      lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem
        (u : 𝒪[K.carrier]) := by
  have hspec := galoisUnitReciprocityMap_spec K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem σ
  rw [hu, unitActionQuotientLift_mk] at hspec
  exact Subtype.ext (Subtype.ext hspec.symm)

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**段 1 のノルム版**:
`σ(α) = [u]_f(α)` かつ `u = 1 + π^i·w`(`w ∈ 𝒪_K^×`)ならば

    `‖σ·α − α‖ = ‖α‖^{q^i}`.

1c(`norm_lubinTateAction_one_add_sub`)で `‖σ·α − α‖ = ‖[π^i w](α)‖` にし、
1a(`action_pi_pow_mul_unit_mem_iteratedLubinTatePsiTorsionPoints`)で
`[π^i w](α)` が `ψ_{k+1}` の根であることを出し、
`spectralNorm_of_mem_iteratedLubinTatePsiTorsionPoints` で両者のノルムを
`‖π‖` の冪で書いて、抽象核 `pow_rpow_one_div_mul` + `nat_pow_sub_factor` で合わせる。 -/
theorem norm_smul_torsionGen_sub
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (i k : ℕ) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1)
      (by omega))
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1))
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (σ : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≃ₐ[K.carrier]
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (u : (𝒪[K.carrier])ˣ) (w : (𝒪[K.carrier])ˣ)
    (hval : (u : 𝒪[K.carrier]) = 1 + π ^ i * (w : 𝒪[K.carrier]))
    (hσ : σ • torsionGen K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem =
      lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem
        (u : 𝒪[K.carrier])) :
    ‖σ • torsionGen K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem -
        torsionGen K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem‖ =
      ‖torsionGen K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem‖ ^ ((pp ^ ff) ^ i) := by
  have hq2 : 2 ≤ pp ^ ff := by rw [← hq]; exact Fintype.one_lt_card
  have hone : lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem 1 =
      torsionGen K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem :=
    lubinTateActionAtTorsionPoint_one K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem
  -- 1c: `‖σ·α − α‖ = ‖[π^i w](α)‖`
  have h1c : ‖σ • torsionGen K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem -
      torsionGen K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem‖ =
      ‖lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem
        (π ^ i * (w : 𝒪[K.carrier]))‖ := by
    rw [hσ, ← hone, hval]
    exact norm_lubinTateAction_one_add_sub K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem
      (π ^ i * (w : 𝒪[K.carrier]))
  -- 1a: `[π^i w](α)` は `ψ_{k+1}` の根
  have h1a := action_pi_pow_mul_unit_mem_iteratedLubinTatePsiTorsionPoints
    K hq hπmax hπne0 f hf0 hf1 hf i k x hxψ hxn hmem w
  have hβ : ‖lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem
      (π ^ i * (w : 𝒪[K.carrier]))‖ =
      ‖π‖ ^ (1 / (((pp ^ ff) ^ (k + 1) - (pp ^ ff) ^ k : ℕ) : ℝ)) := by
    show spectralNorm K.carrier K.closure
      (↑(↑(lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem
          (π ^ i * (w : 𝒪[K.carrier])) :
          IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)) = _
    exact spectralNorm_of_mem_iteratedLubinTatePsiTorsionPoints
      K hq hπmax hπne0 f hf0 hf1 hf (k + 1) (by omega) _ h1a
  have hα : ‖torsionGen K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem‖ =
      ‖π‖ ^ (1 / (((pp ^ ff) ^ (i + k + 1) - (pp ^ ff) ^ (i + k) : ℕ) : ℝ)) := by
    show spectralNorm K.carrier K.closure
      (↑(↑(torsionGen K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem :
          IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)) = _
    exact spectralNorm_of_mem_iteratedLubinTatePsiTorsionPoints
      K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) (by omega) x hxψ
  have hdne : ((pp ^ ff) ^ (k + 1) - (pp ^ ff) ^ k : ℕ) ≠ 0 := by
    have : (pp ^ ff) ^ k < (pp ^ ff) ^ (k + 1) :=
      Nat.pow_lt_pow_right (by omega) (by omega)
    omega
  have hNne : ((pp ^ ff) ^ i : ℕ) ≠ 0 := by positivity
  rw [h1c, hβ, hα, nat_pow_sub_factor (pp ^ ff) i k]
  exact (pow_rpow_one_div_mul (norm_nonneg (π : K.carrier)) ((pp ^ ff) ^ i)
    ((pp ^ ff) ^ (k + 1) - (pp ^ ff) ^ k) hNne hdne).symm

/-- **`v_K(u−1) = i` の主単数語での言い換え**: `u ∈ 1+𝔭^i` かつ `u ∉ 1+𝔭^{i+1}` ならば
`u = 1 + π^i·w` と単数 `w` で書ける。★消費側(`hρ`)は `principalUnits` の語で書かれて
いるので、この 1 本があると本ファイルの主定理をそのまま当てられる。 -/
theorem exists_unit_of_mem_principalUnits_sdiff
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (i : ℕ) (u : (𝒪[K.carrier])ˣ)
    (h : u ∈ principalUnits K π i) (hnot : u ∉ principalUnits K π (i + 1)) :
    ∃ w : (𝒪[K.carrier])ˣ, (u : 𝒪[K.carrier]) = 1 + π ^ i * (w : 𝒪[K.carrier]) := by
  obtain ⟨c, hc⟩ := (mem_principalUnits_iff K π i u).mp h
  have hcu : IsUnit c := by
    by_contra hcu
    have hcmem : c ∈ IsLocalRing.maximalIdeal (𝒪[K.carrier]) := mem_nonunits_iff.mpr hcu
    rw [hπmax, Ideal.mem_span_singleton] at hcmem
    obtain ⟨d, hd⟩ := hcmem
    exact hnot ((mem_principalUnits_iff K π (i + 1) u).mpr ⟨d, by rw [hc, hd]; ring⟩)
  exact ⟨hcu.unit, by rw [hc, IsUnit.unit_spec, mul_comm]⟩

section RamIndex

attribute [local instance] isDiscreteValuationRing_adjoinIntegers

/-- ★★★★★★★★★★★★★★★★**原始捩れ点は `𝒪_{K^m_x}` の素元**(原典の「`α ∈ µ^×_{f,m}`」が
「`K^m_x` の素元」であること)。

木には**完全分岐**(`isTotallyRamifiedAdjoin_iteratedLubinTatePsi`)と**素元の存在**
(`exists_uniformizer_adjoin_eq_adjoinIntegers`)しか無かったので、その 2 つと
ノルムの計算で `x` 自身が素元であることを出す:

* 素元 `ϖ` を 1 つ取ると `v(π) = e = [K^m_x : K] = q^m − q^{m−1}`
  (`addVal_algebraMap_eq_ramificationIndex` + `ramificationIndex_mul_inertiaDegree` +
  `finrank_adjoin_iteratedLubinTatePsi`)、よって `‖π‖ = ‖ϖ‖^e`(橋の `→` 向き)。
* 一方 `‖x‖ = ‖π‖^{1/e}`(`spectralNorm_of_mem_iteratedLubinTatePsiTorsionPoints`)。
* 合わせて `‖x‖ = ‖ϖ‖`、橋の `←` 向きで `v(x) = 1 = v(ϖ)`、
  `addVal_eq_iff_associated` で `x` は `ϖ` に同伴。 -/
theorem irreducible_torsionGen
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (n : ℕ) (hn : 1 ≤ n) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n hn)
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))] :
    Irreducible (torsionGen K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem) := by
  haveI := valuationRing_isDVR K
  have hq2 : 2 ≤ pp ^ ff := by rw [← hq]; exact Fintype.one_lt_card
  set e : ℕ := (pp ^ ff) ^ n - (pp ^ ff) ^ (n - 1) with he_def
  have hene : e ≠ 0 := by
    have hlt : (pp ^ ff) ^ (n - 1) < (pp ^ ff) ^ n := Nat.pow_lt_pow_right (by omega) (by omega)
    omega
  have ht : IsTotallyRamifiedAdjoin K x :=
    isTotallyRamifiedAdjoin_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem
  obtain ⟨ϖ, hϖspan, -⟩ := exists_uniformizer_adjoin_eq_adjoinIntegers K x ht
  have hϖirr : Irreducible ϖ := (IsDiscreteValuationRing.irreducible_iff_uniformizer ϖ).mpr hϖspan
  have hϖpos : 0 < ‖ϖ‖ := norm_pos_iff.mpr hϖirr.ne_zero
  have hπnu : ¬ IsUnit π := by
    have hmemπ : π ∈ IsLocalRing.maximalIdeal (𝒪[K.carrier]) := by
      rw [hπmax]; exact Ideal.mem_span_singleton_self π
    exact mem_nonunits_iff.mp hmemπ
  -- `v(π) = e`
  have hπirr : Irreducible π := (IsDiscreteValuationRing.irreducible_iff_uniformizer π).mpr hπmax
  have hram : ramificationIndex K x = e := by
    have h1 := ramificationIndex_mul_inertiaDegree K x
    rw [show inertiaDegree K x = 1 from ht, mul_one] at h1
    rw [h1, finrank_adjoin_iteratedLubinTatePsi K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ hxn hmem]
  have hvπ : IsDiscreteValuationRing.addVal (adjoinIntegers K x)
      (algebraMap (𝒪[K.carrier]) (adjoinIntegers K x) π) = (e : ℕ∞) := by
    rw [addVal_algebraMap_eq_ramificationIndex K x hπirr, hram]
  have hNπ : ‖algebraMap (𝒪[K.carrier]) (adjoinIntegers K x) π‖ = ‖(π : K.carrier)‖ := by
    show spectralNorm K.carrier K.closure (algebraMap K.carrier K.closure (π : K.carrier)) = _
    exact spectralNorm_extends _
  have hπnorm : ‖(π : K.carrier)‖ = ‖ϖ‖ ^ e := by
    rw [← hNπ]
    exact norm_eq_pow_of_addVal_eq (N := fun y : adjoinIntegers K x => ‖y‖)
      (fun a => norm_nonneg a) (fun a b => norm_mul a b) (fun a => a.2) hϖirr hϖpos hvπ
  have hπlt : ‖(π : K.carrier)‖ < 1 := by
    rcases lt_or_eq_of_le (Valued.integer.norm_le_one π) with h | h
    · exact h
    · exact absurd (Valued.integer.isUnit_iff_norm_eq_one.mpr h) hπnu
  have hϖlt : ‖ϖ‖ < 1 := by
    have hϖle : ‖ϖ‖ ≤ 1 := ϖ.2
    rcases lt_or_eq_of_le hϖle with h | h
    · exact h
    · rw [h, one_pow] at hπnorm
      rw [hπnorm] at hπlt
      exact absurd hπlt (lt_irrefl 1)
  -- `‖x‖ = ‖ϖ‖`
  have hxnorm : ‖torsionGen K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem‖ = ‖ϖ‖ := by
    have hval : spectralNorm K.carrier K.closure x = ‖(π : K.carrier)‖ ^ (1 / (e : ℝ)) :=
      spectralNorm_of_mem_iteratedLubinTatePsiTorsionPoints
        K hq hπmax hπne0 f hf0 hf1 hf n hn x hxψ
    show spectralNorm K.carrier K.closure
      (↑(↑(torsionGen K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem :
          IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)) = _
    rw [show (↑(↑(torsionGen K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)) = x from rfl,
      hval, hπnorm, rpow_one_div_natCast_pow (norm_nonneg ϖ) e hene]
  have hv1 : IsDiscreteValuationRing.addVal (adjoinIntegers K x)
      (torsionGen K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem) = ((1 : ℕ) : ℕ∞) :=
    addVal_eq_of_norm_eq_pow (N := fun y : adjoinIntegers K x => ‖y‖)
      norm_zero (fun a => norm_nonneg a) (fun a b => norm_mul a b) (fun a => a.2)
      hϖirr hϖpos hϖlt (by rw [hxnorm, pow_one])
  have hassoc : Associated (torsionGen K hq hπmax hπne0 f hf0 hf1 hf n x hxn hmem) ϖ := by
    rw [← IsDiscreteValuationRing.addVal_eq_iff_associated, hv1,
      IsDiscreteValuationRing.addVal_uniformizer hϖirr]
    rfl
  exact hassoc.symm.irreducible hϖirr

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**Yoshida 2008 Proposition 6.14 の段 1**:

    `i(σ) = v(σ(α) − α) = q^{v_K(u−1)}`

`σ(α) = [u]_f(α)`(原典の `ρ_{f,m}(σ) = u`)で `u − 1 = π^i·w`(`w ∈ 𝒪_K^×`)、
`m = i+k+1` のとき、下付き分岐の指数 `i(σ) = ramIndex α σ` はちょうど `q^i`。

★**等号である。**★`α` は捩れ点 `x` そのもの(`torsionGen`)であり、それが
`adjoinIntegers K x` の素元であることは `irreducible_torsionGen` で示してある
(仮定ではない)。 -/
theorem ramIndex_torsionGen_eq_pow
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (i k : ℕ) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1)
      (by omega))
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1))
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (σ : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≃ₐ[K.carrier]
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (u w : (𝒪[K.carrier])ˣ)
    (hval : (u : 𝒪[K.carrier]) = 1 + π ^ i * (w : 𝒪[K.carrier]))
    (hσ : σ • torsionGen K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem =
      lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem
        (u : 𝒪[K.carrier])) :
    ramIndex (torsionGen K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem) σ =
      (((pp ^ ff) ^ i : ℕ) : ℕ∞) := by
  set α := torsionGen K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem with hα_def
  have hq2 : 2 ≤ pp ^ ff := by rw [← hq]; exact Fintype.one_lt_card
  have hnorm := norm_smul_torsionGen_sub K hq hπmax hπne0 f hf0 hf1 hf i k x hxψ hxn hmem σ u w
    hval hσ
  have hαpos : 0 < ‖α‖ := by
    have hval' := spectralNorm_of_mem_iteratedLubinTatePsiTorsionPoints
      K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) (by omega) x hxψ
    have hαeq : ‖α‖ = ‖π‖ ^ (1 / (((pp ^ ff) ^ (i + k + 1) - (pp ^ ff) ^ (i + k) : ℕ) : ℝ)) := by
      show spectralNorm K.carrier K.closure
        (↑(↑(α : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)) = _
      exact hval'
    rw [hαeq]
    have hπpos : (0 : ℝ) < ‖(π : K.carrier)‖ := by
      rw [norm_pos_iff]
      exact fun h => hπne0 (Subtype.ext h)
    exact Real.rpow_pos_of_pos hπpos _
  have hαlt : ‖α‖ < 1 := by
    show spectralNorm K.carrier K.closure
      (↑(↑(α : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) : K.closure)) < 1
    exact spectralNorm_lt_one_of_mem_iteratedLubinTateTorsionPoints
      K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn
  have hirr : Irreducible α :=
    irreducible_torsionGen K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) (by omega) x hxψ hxn hmem
  show IsDiscreteValuationRing.addVal (adjoinIntegers K x) (σ • α - α) = _
  exact addVal_eq_of_norm_eq_pow (N := fun y : adjoinIntegers K x => ‖y‖)
    norm_zero (fun a => norm_nonneg a) (fun a b => norm_mul a b) (fun a => a.2)
    hirr hαpos hαlt hnorm

def ramIndex_torsionGen_eq_pow.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Proposition 6.14", sectionId := "prop-6-14" }

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**段 1(消費側の語で書いた形)**:
`ρ_{f,m}(σ) = u mod 𝔭^m` で `u ∈ 1+𝔭^i \ (1+𝔭^{i+1})`(＝ `v_K(u−1) = i`)ならば

    `i(σ) = q^i`.

`Found/PGC/LubinTateUpperRamificationVanish.lean` の `hρ` は `principalUnits` の語で
書かれているので、`hρ` を供給する側はこの形を使うのがよい。 -/
theorem ramIndex_torsionGen_eq_pow_of_principalUnits
    {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
    (i k : ℕ) (x : K.closure)
    (hxψ : x ∈ iteratedLubinTatePsiTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1)
      (by omega))
    (hxn : x ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1))
    (hmem : x ∈ IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (σ : IntermediateField.adjoin K.carrier ({x} : Set K.closure) ≃ₐ[K.carrier]
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))
    (u : (𝒪[K.carrier])ˣ)
    (hu : galoisUnitReciprocityMap K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) (by omega) x
        hxψ hxn hmem σ = (u : (𝒪[K.carrier])ˣ ⧸ principalUnits K π (i + k + 1)))
    (hmemi : u ∈ principalUnits K π i) (hnot : u ∉ principalUnits K π (i + 1)) :
    ramIndex (torsionGen K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) x hxn hmem) σ =
      (((pp ^ ff) ^ i : ℕ) : ℕ∞) := by
  obtain ⟨w, hw⟩ := exists_unit_of_mem_principalUnits_sdiff K hπmax i u hmemi hnot
  exact ramIndex_torsionGen_eq_pow K hq hπmax hπne0 f hf0 hf1 hf i k x hxψ hxn hmem σ u w hw
    (smul_torsionGen_eq_lubinTateAction K hq hπmax hπne0 f hf0 hf1 hf (i + k + 1) (by omega) x
      hxψ hxn hmem σ u hu)

end RamIndex

end ABC3.Found.PGC
