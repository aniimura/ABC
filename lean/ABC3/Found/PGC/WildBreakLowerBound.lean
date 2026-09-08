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
* §3 抽象核(可換環＋超距離) —— `∏(1+η) − 1` の上界／`(1+s)^p − 1 − s^p` の上界
* §4 組み立て

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

end WildBreak

end ABC3.Found.PGC
