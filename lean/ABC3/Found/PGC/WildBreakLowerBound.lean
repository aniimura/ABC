import ABC3.Found.PGC.TotallyRamifiedLayer

/-!
# [pGC] 暴分岐の跳びの**下界** `1 ≤ u` —— 剰余体を作らずにノルムだけで出す

`Found/PGC/TotallyRamifiedLayer.lean` の出口
`..._of_uniformizer_deg_p` に残った★唯一の分岐の仮説★

    hbreak : ‖g π − π‖ = ‖π‖^{t+1} → 1 ≤ t ∧ (p−1)·t ≤ p·e

のうち、★**下界 `1 ≤ t`**(＝ 次数 `p` の全分岐拡大は**暴**分岐である)を埋める。

## ★★前波(＝自分)の記述の訂正

前波の報告と `TotallyRamifiedLayer.lean` の冒頭 docstring はこう書いた:

> ★下界 `1 ≤ u` の証明には**剰余体**(`ū^p = 1 ⇒ ū = 1`)と `g` が等長であることが要る。
> ★本ファイルの設定(`NormedField M` と `g : M ≃ₐ[K] M`)は `‖g z‖ = ‖z‖` を**含まない**ので、そこで止まる。

★**半分は誤りだった**(★どこが誤りかを名指しする。★その docstring は他が読んでいるので直さない):

1. ★**「剰余体が要る」は誤り。** §2 `exists_sub_algebraMap_norm_le` が示す通り、
   「`‖z‖ ≤ 1` なら `K` の元で `‖π‖` 以内に近似できる」は
   ★**素元の冪が基底であること 1 本**から出る。剰余体は一度も作らない。
2. ★**「等長性を含まない」は正しい**(仮説として要る)。★ただし**無料で供給される**:
   `PureStepSetup.norm_algEquiv_eq`(実測: `lean/ABC3/Found/PGC/PureStepSetup.lean:283`)が
   「`NormedAlgebra k M` のノルムが完備な `k` 上のスペクトルノルムなら `k`-代数同型は等長」を
   与える。★本ファイルは抽象のままなので `hiso` を**仮説で受ける**が、
   `PAdicLocalField` に代入する側では `norm_algHom_eq_of_pAdicLocalField` で埋まる。
   ⇒ ★「等長性で止まる」ではなく「等長性は仮説 1 本で、供給は在庫に在る」が正しい。

## 証明の筋(★圧縮で消えても復元できるようここに書く)

`u := g π / π` と置く。`hiso` から `‖u‖ = 1`。`g^p = 1` の telescoping で

    ∏_{j<p} g^j(u) = ∏_{j<p} (g^{j+1}π / g^jπ) = g^pπ / π = 1.

§2 で `a ∈ K` を取り `‖u − a‖ ≤ ‖π‖`。`‖a‖ = ‖u‖ = 1`。`g^j a = a` なので
`‖g^j u − a‖ = ‖g^j(u − a)‖ = ‖u − a‖ ≤ ‖π‖`。ゆえに `g^j u = a(1 + η_j)`,
`‖η_j‖ ≤ ‖π‖`。§3 の積の核から `‖∏(1+η_j) − 1‖ ≤ ‖π‖`、よって `‖a^p − 1‖ ≤ ‖π‖`。
`a^p − 1 ∈ K` は値群が `‖π‖^{pℤ}` なので `‖a^p − 1‖ ≤ ‖π‖^p`。
§3 の二項の核 `‖(1+s)^p − 1 − s^p‖ ≤ ‖p‖`(`s := a − 1`)と `‖p‖ = ‖π‖^{pe} ≤ ‖π‖^p` から
`‖a−1‖^p = ‖(a−1)^p‖ ≤ ‖π‖^p`、すなわち `‖a − 1‖ ≤ ‖π‖`。
最後に `‖u − 1‖ ≤ max(‖u−a‖, ‖a−1‖) ≤ ‖π‖` で `‖gπ − π‖ = ‖π‖·‖u−1‖ ≤ ‖π‖^2`。★これが `1 ≤ t`。

## 節の構成

* §1 抽象核(超距離だけ) —— 有限和のノルム上界／狭義上界／★**ノルムが相異なる族では各項 ≤ 和**
* §2 ★★★`exists_sub_algebraMap_norm_le` —— **全分岐 ⇒ 剰余体が下りる**(剰余体を作らない)
* §3 抽象核(体＋超距離) —— `∏(1+η) − 1` の上界／`(1+s)^p − 1 − s^p` の上界／
  ★`norm_sub_one_pow_le`(「標数 `p` で `x^p = 1 ⇒ x = 1`」のノルム版)
* §4 ★★★★★★`norm_sub_le_sq_of_totallyRamified` —— **`‖gπ − π‖ ≤ ‖π‖²`**(＝ `1 ≤ t`)
* §5 ★★★★★`..._of_uniformizer_deg_p_upper` —— 出口の分岐の仮説が★**上界 1 本だけ**になった

## ★★何が閉じ、何が残ったか(実測)

| 出口 | 分岐についての仮説 |
|---|---|
| `TotallyRamifiedLayer.…_of_uniformizer_deg_p`(前波) | `1 ≤ t` **かつ** `(p−1)t ≤ p·e` |
| ★本ファイル `…_of_uniformizer_deg_p_upper` | ★**`(p−1)t ≤ p·e` だけ**(＋ `hiso`) |

★残り 1 点は**上界** `(p−1)t ≤ p·e`(Serre, Local Fields IV の「異なり」の評価)。
★その筋は「`f'(π) = ∏_{j≠0}(π − g^jπ)` と Eisenstein 多項式の係数評価」で、
★下界と違って**最小多項式の係数が整**であることが要る。★本ファイルでは**測っていない**。

## ★在庫の測定(コマンドを残す)

```
grep -n "IsUltrametricDist" .cache/mathlib-index.txt | grep -i "sum\|max\|add"
  → norm_add_le_max / norm_add_eq_max_of_norm_ne_norm は在る(★索引には乗法形しか出ない。
    to_additive 生成名なので「無い」は嘘。#既知の 7 例に追加)。
  → ★IsUltrametricDist.norm_sum_le は**無い**(#check で確認)。§1 で自作した。
#check @IsUltrametricDist.norm_natCast_le_one  → 在る (R n : ‖↑n‖ ≤ 1)
#check @zpow_le_zpow_right_of_le_one₀ / @one_lt_zpow_of_neg₀ / @zpow_right_inj₀  → 3 本とも在る
#check @Submodule.mem_span_range_iff_exists_fun → 在る(基底展開はこれ 1 本)
#check @Finset.exists_max_image / @Finset.add_sum_erase → 在る
grep -n "pow_lt_pow_left" .cache/mathlib-index.txt
  → ★`pow_lt_pow_left` は Unknown identifier。正しい名は `pow_lt_pow_left₀` と
    ★`le_of_pow_le_pow_left₀ (hn : n ≠ 0) (hb : 0 ≤ b) : a^n ≤ b^n → a ≤ b`(こちらを使った)。
#check @AlgEquiv.mul_apply       → 在る。`(e₁ * e₂) x = e₁ (e₂ x)`。`g^{j+1}` の分解に使う。
#check @pow_le_one₀ / @pow_le_pow_of_le_one / @zpow_lt_zpow_right_of_lt_one₀ → 3 本とも在る
#check @Nat.Prime.dvd_choose_self → 在る (k ≠ 0 → k < p → p ∣ p.choose k)
```
-/

namespace ABC3.Found.PGC

namespace WildBreak

open Finset

/-! ## §1 抽象核(超距離だけ) -/

section Ultra

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- 各項が `r` 以下なら有限和も `r` 以下(超距離)。 -/
theorem norm_sum_le_of_forall_le {ι : Type*} (s : Finset ι) (x : ι → M) {r : ℝ}
    (hr : 0 ≤ r) (h : ∀ i ∈ s, ‖x i‖ ≤ r) : ‖∑ i ∈ s, x i‖ ≤ r := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using hr
  | insert a s ha ih =>
      rw [Finset.sum_insert ha]
      refine (IsUltrametricDist.norm_add_le_max _ _).trans ?_
      exact max_le (h a (Finset.mem_insert_self a s))
        (ih (fun i hi => h i (Finset.mem_insert_of_mem hi)))

/-- 各項が `r` 未満なら有限和も `r` 未満(超距離、`0 < r`)。 -/
theorem norm_sum_lt_of_forall_lt {ι : Type*} (s : Finset ι) (x : ι → M) {r : ℝ}
    (hr : 0 < r) (h : ∀ i ∈ s, ‖x i‖ < r) : ‖∑ i ∈ s, x i‖ < r := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using hr
  | insert a s ha ih =>
      rw [Finset.sum_insert ha]
      refine lt_of_le_of_lt (IsUltrametricDist.norm_add_le_max _ _) ?_
      exact max_lt (h a (Finset.mem_insert_self a s))
        (ih (fun i hi => h i (Finset.mem_insert_of_mem hi)))

/-- ★★**抽象核** —— ノルムが**相異なる**族では、各項のノルムは**和のノルム以下**。

★超距離で「打ち消しが起きない」ことの中身。分岐・付値の語彙が 1 語も出ない。 -/
theorem norm_le_norm_sum_of_pairwise_ne {ι : Type*} [Fintype ι] [DecidableEq ι] (x : ι → M)
    (hpair : ∀ i j, i ≠ j → ‖x i‖ ≠ ‖x j‖) (i₀ : ι) : ‖x i₀‖ ≤ ‖∑ i, x i‖ := by
  obtain ⟨i₁, -, hmax⟩ :=
    Finset.exists_max_image (Finset.univ : Finset ι) (fun i => ‖x i‖) ⟨i₀, Finset.mem_univ i₀⟩
  rcases (norm_nonneg (x i₁)).lt_or_eq with hpos | hzero
  · have hlt : ‖∑ i ∈ Finset.univ.erase i₁, x i‖ < ‖x i₁‖ := by
      refine norm_sum_lt_of_forall_lt _ _ hpos (fun i hi => ?_)
      have hne : i ≠ i₁ := Finset.ne_of_mem_erase hi
      exact lt_of_le_of_ne (hmax i (Finset.mem_univ i)) (hpair i i₁ hne)
    have hsum : ‖∑ i, x i‖ = ‖x i₁‖ := by
      rw [← Finset.add_sum_erase _ x (Finset.mem_univ i₁),
        IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (ne_of_gt hlt)]
      exact max_eq_left (le_of_lt hlt)
    rw [hsum]
    exact hmax i₀ (Finset.mem_univ i₀)
  · have : ‖x i₀‖ = 0 :=
      le_antisymm (by rw [hzero]; exact hmax i₀ (Finset.mem_univ i₀)) (norm_nonneg _)
    rw [this]
    exact norm_nonneg _

/-- ★★**抽象核(上の一般化)** —— **0 を許した**版: 「0 でない項どうしのノルムが相異なる」だけで
各項のノルムは和のノルム以下。★係数に 0 が混じる基底展開に使う。 -/
theorem norm_le_norm_sum_of_pairwise_ne' {ι : Type*} [Fintype ι] [DecidableEq ι] (x : ι → M)
    (hpair : ∀ i j, i ≠ j → x i ≠ 0 → x j ≠ 0 → ‖x i‖ ≠ ‖x j‖) (i₀ : ι) :
    ‖x i₀‖ ≤ ‖∑ i, x i‖ := by
  obtain ⟨i₁, -, hmax⟩ :=
    Finset.exists_max_image (Finset.univ : Finset ι) (fun i => ‖x i‖) ⟨i₀, Finset.mem_univ i₀⟩
  rcases (norm_nonneg (x i₁)).lt_or_eq with hpos | hzero
  · have hx1 : x i₁ ≠ 0 := norm_pos_iff.mp hpos
    have hlt : ‖∑ i ∈ Finset.univ.erase i₁, x i‖ < ‖x i₁‖ := by
      refine norm_sum_lt_of_forall_lt _ _ hpos (fun i hi => ?_)
      have hne : i ≠ i₁ := Finset.ne_of_mem_erase hi
      rcases eq_or_ne (x i) 0 with h0 | h0
      · rw [h0, norm_zero]; exact hpos
      · exact lt_of_le_of_ne (hmax i (Finset.mem_univ i)) (hpair i i₁ hne h0 hx1)
    have hsum : ‖∑ i, x i‖ = ‖x i₁‖ := by
      rw [← Finset.add_sum_erase _ x (Finset.mem_univ i₁),
        IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (ne_of_gt hlt)]
      exact max_eq_left (le_of_lt hlt)
    rw [hsum]
    exact hmax i₀ (Finset.mem_univ i₀)
  · have : ‖x i₀‖ = 0 :=
      le_antisymm (by rw [hzero]; exact hmax i₀ (Finset.mem_univ i₀)) (norm_nonneg _)
    rw [this]
    exact norm_nonneg _

end Ultra

/-! ## §2 ★★★全分岐なら**剰余体は下りる**(`f = 1`)-/

section Residue

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★★★★**全分岐 ⇒ 剰余体が下りる**。`‖z‖ ≤ 1` なら `a ∈ K` が在って `‖z − a‖ ≤ ‖π‖`。

★証明は「素元の冪が基底」(`TotallyRamified.linearIndependent_of_ne_mod`)と、§1 の
「ノルムが相異なる族では各項 ≤ 和」だけ。★**剰余体そのものを一度も作らない。** -/
theorem exists_sub_algebraMap_norm_le [FiniteDimensional K M] {π : M} {n : ℕ} (hn1 : 1 < n)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    {z : M} (hz : ‖z‖ ≤ 1) : ∃ a : K, ‖z - algebraMap K M a‖ ≤ ‖π‖ := by
  classical
  have hπne : ‖π‖ ≠ 1 := ne_of_lt hπ1
  have hli : LinearIndependent K (fun l : Fin n => π ^ (l : ℕ)) := by
    refine TotallyRamified.linearIndependent_of_ne_mod hπ0 hπne hvalK _
      (fun l => ((l : ℕ) : ℤ)) ?_ ?_
    · intro l; simp [zpow_natCast]
    · intro i j hij
      exact TotallyRamifiedLayer.not_dvd_sub_of_lt i.isLt j.isLt (fun h => hij (Fin.ext h))
  have hzmem : z ∈ Submodule.span K (Set.range (fun l : Fin n => π ^ (l : ℕ))) := by
    rw [hli.span_eq_top_of_card_eq_finrank' (by simp [hn])]; exact Submodule.mem_top
  obtain ⟨c, hc⟩ := (Submodule.mem_span_range_iff_exists_fun K).mp hzmem
  set x : Fin n → M := fun l => c l • π ^ (l : ℕ) with hxdef
  have hxval : ∀ l : Fin n, x l = c l • π ^ (l : ℕ) := fun _ => rfl
  have hsum : ∑ l, x l = z := hc
  have hexp : ∀ l : Fin n, c l ≠ 0 →
      ∃ m : ℤ, ‖x l‖ = ‖π‖ ^ ((n : ℤ) * m + ((l : ℕ) : ℤ)) := by
    intro l hl
    obtain ⟨m, hm⟩ := hvalK (c l) hl
    exact ⟨m, by rw [hxval, Algebra.smul_def, norm_mul, hm, norm_pow, ← zpow_natCast ‖π‖ (l : ℕ),
      ← zpow_add₀ (ne_of_gt hπ0)]⟩
  have hczero : ∀ l : Fin n, x l ≠ 0 → c l ≠ 0 := by
    intro l hx hcl
    exact hx (by rw [hxval, hcl, zero_smul])
  have hpair : ∀ i j : Fin n, i ≠ j → x i ≠ 0 → x j ≠ 0 → ‖x i‖ ≠ ‖x j‖ := by
    intro i j hij hi hj
    obtain ⟨mi, hmi⟩ := hexp i (hczero i hi)
    obtain ⟨mj, hmj⟩ := hexp j (hczero j hj)
    rw [hmi, hmj]
    intro hcon
    have h2 : (n : ℤ) * mi + ((i : ℕ) : ℤ) = (n : ℤ) * mj + ((j : ℕ) : ℤ) :=
      (zpow_right_inj₀ hπ0 hπne).mp hcon
    refine TotallyRamifiedLayer.not_dvd_sub_of_lt i.isLt j.isLt
      (fun h => hij (Fin.ext h)) ⟨mj - mi, ?_⟩
    linarith [mul_sub (n : ℤ) mj mi]
  have hle1 : ∀ l : Fin n, ‖x l‖ ≤ 1 := fun l =>
    le_trans (by rw [← hsum]; exact norm_le_norm_sum_of_pairwise_ne' x hpair l) hz
  obtain ⟨i0, hi0⟩ : ∃ i : Fin n, (i : ℕ) = 0 := ⟨⟨0, by omega⟩, rfl⟩
  refine ⟨c i0, ?_⟩
  have h0 : x i0 = algebraMap K M (c i0) := by
    rw [hxval, Algebra.smul_def, hi0, pow_zero, mul_one]
  have hrest : z - algebraMap K M (c i0) = ∑ l ∈ Finset.univ.erase i0, x l := by
    rw [← h0, ← hsum, ← Finset.add_sum_erase _ x (Finset.mem_univ i0)]
    ring
  rw [hrest]
  refine norm_sum_le_of_forall_le _ _ (le_of_lt hπ0) (fun l hl => ?_)
  rcases eq_or_ne (c l) 0 with hcz | hcz
  · rw [hxval, hcz, zero_smul, norm_zero]; exact le_of_lt hπ0
  · obtain ⟨m, hm⟩ := hexp l hcz
    have hlne : (l : ℕ) ≠ 0 := by
      intro h
      exact (Finset.ne_of_mem_erase hl) (Fin.ext (by rw [h, hi0]))
    have hEpos : 0 ≤ (n : ℤ) * m + ((l : ℕ) : ℤ) := by
      by_contra hcon
      rw [not_le] at hcon
      have hgt : (1:ℝ) < ‖π‖ ^ ((n : ℤ) * m + ((l : ℕ) : ℤ)) :=
        one_lt_zpow_of_neg₀ hπ0 hπ1 hcon
      rw [← hm] at hgt
      linarith [hle1 l]
    have hEne : (n : ℤ) * m + ((l : ℕ) : ℤ) ≠ 0 := by
      intro hE
      refine TotallyRamifiedLayer.not_dvd_sub_of_lt l.isLt (by omega : 0 < n) hlne ⟨-m, ?_⟩
      push_cast
      linarith [mul_neg (n : ℤ) m]
    rw [hm]
    calc ‖π‖ ^ ((n : ℤ) * m + ((l : ℕ) : ℤ)) ≤ ‖π‖ ^ (1 : ℤ) :=
          zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) (by omega)
      _ = ‖π‖ := zpow_one _

end Residue

/-! ## §3 抽象核(体＋超距離) —— 積・二項・「`x^p = 1 ⇒ x = 1`」のノルム版 -/

section ProdBinom

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★★**抽象核** —— `‖η i‖ ≤ r ≤ 1` なら `‖∏ (1 + η i) − 1‖ ≤ r`。 -/
theorem norm_prod_one_add_sub_one_le {ι : Type*} (s : Finset ι) (η : ι → M) {r : ℝ}
    (hr0 : 0 ≤ r) (hr1 : r ≤ 1) (h : ∀ i ∈ s, ‖η i‖ ≤ r) :
    ‖(∏ i ∈ s, (1 + η i)) - 1‖ ≤ r := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using hr0
  | insert a s ha ih =>
      have hha : ‖η a‖ ≤ r := h a (Finset.mem_insert_self a s)
      have h1 : ‖1 + η a‖ ≤ 1 := by
        refine (IsUltrametricDist.norm_add_le_max _ _).trans ?_
        simpa using hha.trans hr1
      have h2 : ‖(∏ i ∈ s, (1 + η i)) - 1‖ ≤ r :=
        ih (fun i hi => h i (Finset.mem_insert_of_mem hi))
      rw [Finset.prod_insert ha,
        show (1 + η a) * (∏ i ∈ s, (1 + η i)) - 1
          = (1 + η a) * ((∏ i ∈ s, (1 + η i)) - 1) + η a by ring]
      refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le ?_ hha)
      rw [norm_mul]
      calc ‖1 + η a‖ * ‖(∏ i ∈ s, (1 + η i)) - 1‖ ≤ 1 * r :=
            mul_le_mul h1 h2 (norm_nonneg _) zero_le_one
        _ = r := one_mul r

/-- ★★**抽象核** —— `‖t‖ ≤ 1` なら `‖(1+t)^p − 1 − t^p‖ ≤ ‖(p : M)‖`。

★中身は「`0 < k < p` なら `p ∣ C(p,k)`」だけ。★超距離で各項を潰す。 -/
theorem norm_one_add_pow_sub_le {p : ℕ} (hp : p.Prime) {t : M} (ht : ‖t‖ ≤ 1) :
    ‖(1 + t) ^ p - 1 - t ^ p‖ ≤ ‖(p : M)‖ := by
  classical
  obtain ⟨q, rfl⟩ : ∃ q, p = q + 1 := ⟨p - 1, by have := hp.pos; omega⟩
  have hbin : (1 + t) ^ (q + 1)
      = ∑ k ∈ Finset.range (q + 1 + 1), t ^ k * ((q + 1).choose k : M) := by
    rw [add_comm (1 : M) t]
    simpa using add_pow t 1 (q + 1)
  have hpeel : (1 + t) ^ (q + 1) - 1 - t ^ (q + 1)
      = ∑ k ∈ Finset.range q, t ^ (k + 1) * ((q + 1).choose (k + 1) : M) := by
    rw [hbin, Finset.sum_range_succ, Finset.sum_range_succ']
    simp
    ring
  rw [hpeel]
  refine norm_sum_le_of_forall_le _ _ (norm_nonneg _) (fun k hk => ?_)
  have hkq : k < q := Finset.mem_range.mp hk
  obtain ⟨d, hd⟩ : (q + 1) ∣ (q + 1).choose (k + 1) :=
    hp.dvd_choose_self (Nat.succ_ne_zero k) (by omega)
  rw [norm_mul, hd]
  push_cast
  rw [norm_mul, norm_pow]
  have h1 : ‖t‖ ^ (k + 1) ≤ 1 := pow_le_one₀ (norm_nonneg _) ht
  have h2 : ‖(d : M)‖ ≤ 1 := IsUltrametricDist.norm_natCast_le_one M d
  calc ‖t‖ ^ (k + 1) * (‖((q : M) + 1)‖ * ‖(d : M)‖) ≤ 1 * (‖((q : M) + 1)‖ * 1) :=
        mul_le_mul h1 (mul_le_mul_of_nonneg_left h2 (norm_nonneg _)) (by positivity) zero_le_one
    _ = ‖((q : M) + 1)‖ := by ring

/-- ★★★**抽象核** —— `‖A‖ = 1`・`‖A^p − 1‖ ≤ c`・`‖(p:M)‖ ≤ c` なら `‖A − 1‖^p ≤ c`。

★これが「標数 `p` で `x^p = 1 ⇒ x = 1`」のノルム版である。★体の分岐も剰余体も出てこない。 -/
theorem norm_sub_one_pow_le {p : ℕ} (hp : p.Prime) {A : M} {c : ℝ}
    (hA : ‖A‖ = 1) (hc : ‖(p : M)‖ ≤ c) (hAp : ‖A ^ p - 1‖ ≤ c) : ‖A - 1‖ ^ p ≤ c := by
  have hs1 : ‖A - 1‖ ≤ 1 := by
    have h := IsUltrametricDist.norm_add_le_max A (-1 : M)
    simpa [sub_eq_add_neg, hA] using h
  have hbin : ‖(1 + (A - 1)) ^ p - 1 - (A - 1) ^ p‖ ≤ ‖(p : M)‖ :=
    norm_one_add_pow_sub_le hp hs1
  rw [show (1 : M) + (A - 1) = A by ring] at hbin
  calc ‖A - 1‖ ^ p = ‖(A - 1) ^ p‖ := (norm_pow _ _).symm
    _ = ‖(A ^ p - 1) + -(A ^ p - 1 - (A - 1) ^ p)‖ := by ring_nf
    _ ≤ max ‖A ^ p - 1‖ ‖-(A ^ p - 1 - (A - 1) ^ p)‖ :=
        IsUltrametricDist.norm_add_le_max _ _
    _ ≤ c := by rw [norm_neg]; exact max_le hAp (hbin.trans hc)

end ProdBinom
/-! ## §4 組み立て -/

section Assembly

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

omit [IsUltrametricDist M] in
/-- 等長性は冪に伝わる。 -/
theorem norm_pow_apply (g : M ≃ₐ[K] M) (hiso : ∀ z : M, ‖g z‖ = ‖z‖) :
    ∀ (j : ℕ) (z : M), ‖(g ^ j) z‖ = ‖z‖ := by
  intro j
  induction j with
  | zero => intro z; simp
  | succ j ih => intro z; rw [pow_succ, AlgEquiv.mul_apply, ih, hiso]

omit [IsUltrametricDist M] in
/-- telescoping —— `∏_{j<n} g^{j+1}π / g^jπ = g^nπ / π`。 -/
theorem prod_telescope (g : M ≃ₐ[K] M) {π : M} (hπ : π ≠ 0)
    (hne : ∀ j : ℕ, (g ^ j) π ≠ 0) (n : ℕ) :
    ∏ j ∈ Finset.range n, ((g ^ (j + 1)) π / (g ^ j) π) = (g ^ n) π / π := by
  induction n with
  | zero => simp [div_self hπ]
  | succ n ih =>
      rw [Finset.prod_range_succ, ih]
      field_simp
      exact mul_div_cancel_left₀ _ (hne n)

/-- ★★★★★★**次数 `p` の全分岐拡大は暴分岐である** —— 跳びの**下界**。

`‖g π − π‖ ≤ ‖π‖^2`、すなわち `‖g π − π‖ = ‖π‖^{t+1}` と書いたときの `1 ≤ t`。

★仮説は「`g^p = 1`」「`g` が等長」「全分岐(`hvalK`)」「`[M:K] = p`」「`‖p‖ = ‖π‖^{p·e}`, `e ≥ 1`」。
★等長性 `hiso` は `PureStepSetup.norm_algEquiv_eq`(スペクトルノルム)が供給する。 -/
theorem norm_sub_le_sq_of_totallyRamified [FiniteDimensional K M] {p e : ℕ} (hp : p.Prime)
    {π : M} (g : M ≃ₐ[K] M) (hgp : g ^ p = 1) (hiso : ∀ z : M, ‖g z‖ = ‖z‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = p)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (he : 0 < e) (heM : ‖(p : M)‖ = ‖π‖ ^ (p * e)) :
    ‖g π - π‖ ≤ ‖π‖ ^ 2 := by
  classical
  have hπne : π ≠ 0 := norm_pos_iff.mp hπ0
  have hgj : ∀ (j : ℕ) (z : M), ‖(g ^ j) z‖ = ‖z‖ := norm_pow_apply g hiso
  have hfne : ∀ j : ℕ, (g ^ j) π ≠ 0 := fun j => norm_pos_iff.mp (by rw [hgj]; exact hπ0)
  set u : M := g π / π with hudef
  have hu1 : ‖u‖ = 1 := by rw [hudef, norm_div, hiso]; exact div_self (ne_of_gt hπ0)
  have hgu : ∀ j : ℕ, (g ^ j) u = (g ^ (j + 1)) π / (g ^ j) π := by
    intro j
    rw [hudef, map_div₀, ← AlgEquiv.mul_apply, ← pow_succ]
  have hprod : ∏ j ∈ Finset.range p, (g ^ j) u = 1 := by
    simp_rw [hgu]
    rw [prod_telescope g hπne hfne p, hgp]
    simpa using div_self hπne
  obtain ⟨a, ha⟩ := exists_sub_algebraMap_norm_le (n := p) hp.one_lt hπ0 hπ1 hn hvalK
    (le_of_eq hu1)
  set A : M := algebraMap K M a with hAdef
  have hAnorm : ‖A‖ = 1 := by
    have hlt : ‖A - u‖ < ‖u‖ := by
      rw [norm_sub_rev, hu1]
      exact lt_of_le_of_lt ha hπ1
    have := IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (x := u) (y := A - u)
      (ne_of_gt hlt)
    rw [show u + (A - u) = A by ring] at this
    rw [this, max_eq_left (le_of_lt hlt), hu1]
  have hAne : A ≠ 0 := norm_pos_iff.mp (by rw [hAnorm]; norm_num)
  have hgA : ∀ j : ℕ, (g ^ j) A = A := fun j => (g ^ j).commutes a
  set η : ℕ → M := fun j => (g ^ j) u / A - 1 with hηdef
  have hfactor : ∀ j : ℕ, (g ^ j) u = A * (1 + η j) := by
    intro j
    rw [hηdef]
    field_simp
    ring
  have hηnorm : ∀ j : ℕ, ‖η j‖ ≤ ‖π‖ := by
    intro j
    have h1 : η j = ((g ^ j) u - A) / A := by rw [hηdef]; field_simp
    have h2 : (g ^ j) u - A = (g ^ j) (u - A) := by rw [map_sub, hgA j]
    rw [h1, norm_div, hAnorm, div_one, h2, hgj]
    exact ha
  have hprodA : A ^ p * ∏ j ∈ Finset.range p, (1 + η j) = 1 := by
    have hcong : ∏ j ∈ Finset.range p, (g ^ j) u = ∏ j ∈ Finset.range p, (A * (1 + η j)) :=
      Finset.prod_congr rfl (fun j _ => hfactor j)
    rw [Finset.prod_mul_distrib, Finset.prod_const, Finset.card_range] at hcong
    rw [← hcong]
    exact hprod
  have hP : ‖(∏ j ∈ Finset.range p, (1 + η j)) - 1‖ ≤ ‖π‖ :=
    norm_prod_one_add_sub_one_le _ _ (le_of_lt hπ0) (le_of_lt hπ1) (fun j _ => hηnorm j)
  have hApm1 : ‖A ^ p - 1‖ ≤ ‖π‖ := by
    have hAp : ‖A ^ p‖ = 1 := by rw [norm_pow, hAnorm, one_pow]
    have hid : A ^ p - 1 = A ^ p * (1 - ∏ j ∈ Finset.range p, (1 + η j)) := by
      rw [mul_sub, mul_one, hprodA]
    rw [hid, norm_mul, hAp, one_mul, ← norm_neg]
    simpa using hP
  -- `A^p − 1` は `K` の元なので、値群から `‖π‖^p` まで下がる
  have hApK : A ^ p - 1 = algebraMap K M (a ^ p - 1) := by
    rw [map_sub, map_pow, map_one, hAdef]
  have hApm1' : ‖A ^ p - 1‖ ≤ ‖π‖ ^ p := by
    rcases eq_or_ne (a ^ p - 1) 0 with h0 | h0
    · rw [hApK, h0, map_zero, norm_zero]
      exact pow_nonneg (norm_nonneg π) p
    · obtain ⟨m, hm⟩ := hvalK (a ^ p - 1) h0
      rw [hApK, hm]
      have hlt1 : ‖π‖ ^ ((p : ℤ) * m) < 1 := by
        rw [← hm, ← hApK]
        exact lt_of_le_of_lt hApm1 hπ1
      have hpos : 0 < (p : ℤ) * m := (zpow_lt_one_iff_right_of_lt_one₀ hπ0 hπ1).mp hlt1
      have hm1 : 1 ≤ m := by
        by_contra hcon
        rw [not_le] at hcon
        nlinarith [hpos, (by exact_mod_cast hp.pos : (0:ℤ) < p)]
      have hple : (p : ℤ) ≤ (p : ℤ) * m := by
        nlinarith [(by exact_mod_cast hp.pos : (0:ℤ) < p)]
      calc ‖π‖ ^ ((p : ℤ) * m) ≤ ‖π‖ ^ ((p : ℕ) : ℤ) :=
            zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) hple
        _ = ‖π‖ ^ p := zpow_natCast _ _
  -- `‖p‖ ≤ ‖π‖^p`
  have hcp : ‖(p : M)‖ ≤ ‖π‖ ^ p := by
    rw [heM]
    exact pow_le_pow_of_le_one (norm_nonneg π) (le_of_lt hπ1) (Nat.le_mul_of_pos_right p he)
  -- 核: `‖A − 1‖^p ≤ ‖π‖^p`
  have hAB : ‖A - 1‖ ^ p ≤ ‖π‖ ^ p := norm_sub_one_pow_le hp hAnorm hcp hApm1'
  have hA1 : ‖A - 1‖ ≤ ‖π‖ := le_of_pow_le_pow_left₀ hp.pos.ne' (norm_nonneg π) hAB
  -- `‖u − 1‖ ≤ ‖π‖`
  have hu1' : ‖u - 1‖ ≤ ‖π‖ := by
    have hsplit : u - 1 = (u - A) + (A - 1) := by ring
    rw [hsplit]
    exact (IsUltrametricDist.norm_add_le_max _ _).trans (max_le ha hA1)
  -- 仕上げ
  have hfin : g π - π = π * (u - 1) := by
    rw [hudef]
    field_simp
  rw [hfin, norm_mul, sq]
  exact mul_le_mul_of_nonneg_left hu1' (norm_nonneg π)


end Assembly

/-! ## §5 ★★★出口 —— 残る分岐の仮説は**上界だけ**になった -/

section Exit

open GainedTowerModel GainedTowerModel.GaloisTower TotallyRamifiedLayer

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★★★★★**`..._of_uniformizer_deg_p` の `hbreak` から下界が消えた形**。

`TotallyRamifiedLayer.exists_norm_sub_algebraMap_le_axDecay_of_uniformizer_deg_p` は
分岐の仮説として `1 ≤ t ∧ (p−1)·t ≤ p·e` を要求していた。
★本定理は**下界 `1 ≤ t` を証明して消し**、残りを `hupper`(上界)1 本にする。
★代償は `hiso`(`g` が等長)で、これは `PureStepSetup.norm_algEquiv_eq` が供給する。 -/
theorem exists_norm_sub_algebraMap_le_axDecay_of_uniformizer_deg_p_upper
    {K : Type*} [Field K] [Algebra K M] [FiniteDimensional K M]
    {p : ℕ} [Fact p.Prime] {π : M}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p) (hiso : ∀ z : M, ‖g z‖ = ‖z‖)
    (hupper : ∀ (e : ℕ) (t : ℤ), 0 < e → ‖(p : M)‖ = ‖π‖ ^ (p * e) →
      ‖g π - π‖ = ‖π‖ ^ (t + 1) → ((p : ℤ) - 1) * t ≤ (p : ℤ) * (e : ℤ))
    (hnK : Module.finrank K M = p)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹) (hπlt : ‖π‖ < 1) (x : M) :
    ∃ y : twr g p 0, ‖x - algebraMap (twr g p 0) M y‖ ≤ axDecay p 1 * ‖g x - x‖ := by
  have hp : p.Prime := Fact.out
  have hπ0 : 0 < ‖π‖ := norm_pos_of_valK (p := p) (n := p) hnormp hvalK
  refine exists_norm_sub_algebraMap_le_axDecay_of_uniformizer_deg_p g hg ?_ hnK hvalK hnormp hπlt x
  intro e t he heM hbr
  refine ⟨?_, hupper e t he heM hbr⟩
  have hgp : g ^ p = 1 := by rw [← hg]; exact pow_orderOf_eq_one g
  have hle := norm_sub_le_sq_of_totallyRamified hp g hgp hiso hπ0 hπlt hnK hvalK he heM
  by_contra hcon
  rw [not_le] at hcon
  have hz2 : ‖π‖ ^ (2 : ℤ) = ‖π‖ ^ (2 : ℕ) := by
    rw [show (2 : ℤ) = ((2 : ℕ) : ℤ) by norm_num, zpow_natCast]
  have hlt2 : ‖π‖ ^ (2 : ℤ) < ‖π‖ ^ (t + 1) :=
    zpow_lt_zpow_right_of_lt_one₀ hπ0 hπlt (by omega)
  rw [hz2, ← hbr] at hlt2
  linarith

end Exit

/-! ## §6 `.src`(原典の対応箇所) -/

def exists_sub_algebraMap_norm_le.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_sub_le_sq_of_totallyRamified.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_axDecay_of_uniformizer_deg_p_upper.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §7 使っている公理の一覧 -/

#print axioms norm_sum_le_of_forall_le
#print axioms norm_le_norm_sum_of_pairwise_ne
#print axioms norm_le_norm_sum_of_pairwise_ne'
#print axioms exists_sub_algebraMap_norm_le
#print axioms norm_prod_one_add_sub_one_le
#print axioms norm_one_add_pow_sub_le
#print axioms norm_sub_one_pow_le
#print axioms norm_sub_le_sq_of_totallyRamified
#print axioms exists_norm_sub_algebraMap_le_axDecay_of_uniformizer_deg_p_upper

end WildBreak

end ABC3.Found.PGC
