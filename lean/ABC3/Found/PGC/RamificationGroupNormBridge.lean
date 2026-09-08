import ABC3.Found.PGC.HasseArfCongruence

/-!
# [pGC] `hrec` の橋 —— **(a) は取れた。止まるのは「環の型を建てる」ところ**

前波(第 1115/1116)で `harith` の (3) Hasse–Arf の合同は
★**残り `hrec` 1 本**まで絞れた:

    hrec : ∀ m, φ (m+1) = φ m + (u (m+1) − u m) / p^{m+1}   (φ m := herbrandPhiGroup G π' (u m))

本ファイルはその橋を (a)(b)(c) に割って測り、★**(a) を閉じ、(c) が在庫であることを確かめ、
(b) が止まる場所を `file:line` で確定した**。

## ★★(a) 単元での 1 つ分の改良 —— **取れた**

前波は `JumpMono.norm_sub_apply_le_mul`(`‖h z − z‖ ≤ ‖z‖·‖π‖^t`)しか持たず、
`‖z‖ ≤ 1` では `‖π‖^t` で ★**`h ∈ G_{t−1}` までしか出なかった**(真は `h ∈ G_t`)。

★鍵は `norm_le_norm_pi_of_le_one`: `0 < l < n` で `‖c·π^l‖ ≤ 1` なら `‖c·π^l‖ ≤ ‖π‖`。
指数 `n·m + l` は `≤ 1` から `≥ 0`、`0 < l < n` から `n ∤ l` ゆえ `≠ 0`、よって `≥ 1`。
これで `l ≥ 1` の項が `‖π‖·‖π‖^t = ‖π‖^{t+1}` に収まり、`l = 0` の項は消える。

⇒ ★★`norm_sub_apply_le_of_norm_le_one` と、両向きの
★★★`mem_ramification_iff` : `(∀ z, ‖z‖ ≤ 1 → ‖h z − z‖ ≤ ‖π‖^{i+1}) ↔ i ≤ t`。
★左辺は `LowerRamificationGroup.lean:270 mem_lowerRamificationGroup_iff_forall`
(`σ ∈ G_i ↔ ∀ x : B, σ•x − x ∈ 𝔪^{i+1}`)の**ノルム版**そのものである。

★持ち場が指した `:277 mem_lowerRamificationGroup_iff`(`B = A[α]` の下で `σα − α` だけ見る形)は
**使わなかった**。★上の `:270` の全称形で足り、しかも `hadj` を要求しないからである。

## ★(c) 漸化式 —— **在庫に在った**(測定コマンドつき)

```
sed -n '455,470p' lean/ABC3/Found/PGC/HerbrandComposition.lean
  → herbrandPhiGroup_natCast :
      herbrandPhiGroup G α (n : ℝ) = (∑ i ∈ Finset.Icc 1 n, (Nat.card (lowerRamificationGroup B G i) : ℝ)) / (Nat.card G : ℝ)
sed -n '183,195p' lean/ABC3/Found/PGC/HasseArfStrongInduction.lean
  → herbrandPhi_succ_natCast(1 段版。分子に |H_{n+1}| を足す形)
```
⇒ `φ(n+1) − φ(n) = |G_{n+1}|/|G|` は `Finset.sum_Icc_succ_top` で 1 行。
★**(c) は書く必要が無い。**

## ★★(b) `|G_i|` を `u` で書き下す —— **ここで止まる。止まる理由は数学ではなく型**

`mem_ramification_iff` から、ノルムの言葉では

    g^j ∈ G_i  ⟺  i ≤ (g^j の跳び)  ,  (g^j の跳び) = u_{v_p(j)}

(後半は前波 `JumpMono`/`WildBreakUpper.norm_pow_apply_sub_eq` の「生成元が同じなら等距離」)。
よって `|G_i| = #{j < p^{k+1} : v_p(j) ≥ m_i} = p^{k+1−m_i}`(`m_i := min{m : u m ≥ i}`)。
★**数学はこれで尽きている。**

★止まるのは `herbrandPhiGroup` が要求する**型**の側である。代入するには

| 要るもの | 現状(実測) |
|---|---|
| `B` : `CommRing`・`IsDomain`・`IsDiscreteValuationRing` | ★`AdjoinIntegers.lean:70 adjoinIntegers`(`{y : K(x) | ‖y‖ ≤ 1}` の `Subring`)が**型としては在る**が、`IsDiscreteValuationRing` のインスタンスは**同ファイルに無い**(`grep -n IsDiscreteValuationRing AdjoinIntegers.lean` は仮定側の行しか出ない) |
| `MulSemiringAction G B` | ★**無い**。ただし `hiso`(等長)から「`g` は `‖·‖ ≤ 1` を保つ」ので**作れるはず**(本ファイルは作っていない) |
| `maximalIdeal B = Ideal.span {π}` | ★**無い** |
| `Algebra.adjoin 𝒪_K {π} = ⊤`(`𝒪_M = 𝒪_K[π]`) | ★中身は本ファイルの `norm_algebraMap_le_one_of_le_one` + `JumpMono.exists_coeff_norm_le` で**尽きている**(係数が `𝒪_K` に入ることを示した)。型に載せていないだけ |

★★さらに `CLAUDE.md`/`lean-idioms.md #69` が
「`adjoinField` / `adjoinIntegers` の境界は越えられない(212 秒 timeout)」と記録している。
⇒ ★**(b) は「数学の穴」ではなく「`𝒪_M` を型として建てて 4 つのインスタンスを載せる」配管であり、
かつ #69 の危険区間に入る。**★本ファイルはそこに手を付けていない。

## 本ファイルが出したもの

* `norm_le_norm_pi_of_le_one`(指数の議論)
* ★`norm_sub_apply_le_of_norm_le_one`((a) の本体)
* ★★`mem_ramification_iff`(分岐群のノルム版、両向き)
* ★`norm_algebraMap_le_one_of_le_one`(係数も整 ＝ `𝒪_M = 𝒪_K[π]` の中身)

## 逸脱の記録

1. `hiso`(等長)は仮説。`ConcreteDegPFree.hiso_g2` の道(スペクトルノルムなら `rfl` 1 行)で供給される。
2. 既存の `Found/PGC/*.lean` は **1 行も書き換えていない**(import のみ)。
-/

namespace ABC3.Found.PGC

namespace RamNormBridge

open Finset WildBreak JumpMono

section Core

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

omit [IsUltrametricDist M] in
/-- ★**指数の議論** —— `0 < l < n` で `‖c·π^l‖ ≤ 1` なら `‖c·π^l‖ ≤ ‖π‖`。

`‖c·π^l‖ = ‖π‖^{n·m + l}` で、`≤ 1` から指数 `≥ 0`、`0 < l < n` から `n ∤ l` ゆえ指数 `≠ 0`。
★これが「単元の展開では `l = 0` の項だけが大きい」ことの中身である。 -/
theorem norm_le_norm_pi_of_le_one {π : M} {n : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    {c : K} {l : ℕ} (hl0 : 0 < l) (hln : l < n)
    (hle : ‖algebraMap K M c * π ^ l‖ ≤ 1) :
    ‖algebraMap K M c * π ^ l‖ ≤ ‖π‖ := by
  rcases eq_or_ne c 0 with hc0 | hc0
  · rw [hc0, map_zero, zero_mul, norm_zero]; exact le_of_lt hπ0
  · obtain ⟨m, hm⟩ := hvalK c hc0
    have hval : ‖algebraMap K M c * π ^ l‖ = ‖π‖ ^ ((n : ℤ) * m + (l : ℤ)) := by
      rw [norm_mul, hm, norm_pow, ← zpow_natCast ‖π‖ l, ← zpow_add₀ (ne_of_gt hπ0)]
    have hE0 : 0 ≤ (n : ℤ) * m + (l : ℤ) := by
      by_contra hcon
      rw [not_le] at hcon
      have hgt : (1:ℝ) < ‖π‖ ^ ((n : ℤ) * m + (l : ℤ)) := one_lt_zpow_of_neg₀ hπ0 hπ1 hcon
      rw [← hval] at hgt
      linarith
    have hEne : (n : ℤ) * m + (l : ℤ) ≠ 0 := by
      intro hE
      refine TotallyRamifiedLayer.not_dvd_sub_of_lt hln (by omega : 0 < n)
        (by omega : l ≠ 0) ⟨-m, ?_⟩
      push_cast
      linarith [mul_neg (n : ℤ) m]
    rw [hval]
    calc ‖π‖ ^ ((n : ℤ) * m + (l : ℤ)) ≤ ‖π‖ ^ (1 : ℤ) :=
          zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) (by omega)
      _ = ‖π‖ := zpow_one _

/-- ★★★★**(a) 単元での 1 つ分の改良** —— `‖z‖ ≤ 1` なら `‖h z − z‖ ≤ ‖π‖^{t+1}`。

★前波 `JumpMono.norm_sub_apply_le_mul` は `‖z‖·‖π‖^t` までしか出さず、
`‖z‖ ≤ 1` では `‖π‖^t`(1 つ足りない)であった。★`l = 0` の項が消えることと、
`l ≥ 1` の項が `‖π‖` 以下であること(上の指数の議論)で 1 つ取り戻す。

★これで `h ∈ G_t`(真の値)がノルムの言葉で言える。 -/
theorem norm_sub_apply_le_of_norm_le_one [FiniteDimensional K M] {π : M} {n t : ℕ}
    (h : M ≃ₐ[K] M) (hiso : ∀ w : M, ‖h w‖ = ‖w‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hbr : ‖h π - π‖ = ‖π‖ ^ (t + 1)) {z : M} (hz : ‖z‖ ≤ 1) :
    ‖h z - z‖ ≤ ‖π‖ ^ (t + 1) := by
  classical
  obtain ⟨c, hc, hnorm⟩ := exists_coeff_norm_le hπ0 hπ1 hn hvalK z
  have hzz : h z - z = ∑ l : Fin n, (c l • ((h π) ^ (l : ℕ) - π ^ (l : ℕ))) := by
    rw [← hc, map_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl (fun l _ => ?_)
    rw [Algebra.smul_def, Algebra.smul_def, map_mul, AlgEquiv.commutes, map_pow, mul_sub]
  rw [hzz]
  refine norm_sum_le_of_forall_le _ _ (by positivity) (fun l _ => ?_)
  rcases Nat.eq_zero_or_pos (l : ℕ) with hl0 | hl1
  · rw [hl0]
    simp [pow_nonneg (norm_nonneg π) (t + 1)]
  · have hle1 : ‖algebraMap K M (c l) * π ^ (l : ℕ)‖ ≤ 1 := by
      have h0 := hnorm l
      rw [Algebra.smul_def] at h0
      exact le_trans h0 hz
    have hcl : ‖algebraMap K M (c l) * π ^ (l : ℕ)‖ ≤ ‖π‖ :=
      norm_le_norm_pi_of_le_one hπ0 hπ1 hvalK hl1 l.isLt hle1
    have hstep : ‖(h π) ^ (l : ℕ) - π ^ (l : ℕ)‖ ≤ ‖π‖ ^ ((l : ℕ) - 1) * ‖π‖ ^ (t + 1) := by
      rw [← hbr]
      exact GainedTowerStep.norm_pow_sub_pow_le hπ0 (le_of_eq (hiso π)) _
    have hcomb : ‖π‖ ^ ((l : ℕ) - 1) * ‖π‖ ^ (t + 1) = ‖π‖ ^ (l : ℕ) * ‖π‖ ^ t := by
      rw [← pow_add, ← pow_add]
      congr 1
      omega
    rw [Algebra.smul_def, norm_mul]
    calc ‖algebraMap K M (c l)‖ * ‖(h π) ^ (l : ℕ) - π ^ (l : ℕ)‖
        ≤ ‖algebraMap K M (c l)‖ * (‖π‖ ^ (l : ℕ) * ‖π‖ ^ t) := by
          rw [← hcomb]; exact mul_le_mul_of_nonneg_left hstep (norm_nonneg _)
      _ = (‖algebraMap K M (c l)‖ * ‖π‖ ^ (l : ℕ)) * ‖π‖ ^ t := by ring
      _ ≤ ‖π‖ * ‖π‖ ^ t := by
          refine mul_le_mul_of_nonneg_right ?_ (by positivity)
          rwa [← norm_pow, ← norm_mul]
      _ = ‖π‖ ^ (t + 1) := by rw [pow_succ]; ring

/-- ★★★★★**分岐群のノルム言語での特徴づけ**(両向き)。

    (∀ z, ‖z‖ ≤ 1 → ‖h z − z‖ ≤ ‖π‖^{i+1})  ↔  i ≤ t        (`‖hπ − π‖ = ‖π‖^{t+1}`)

★左辺は `LowerRamificationGroup.lean:270` の
`σ ∈ lowerRamificationGroup B G i ↔ ∀ x : B, σ•x − x ∈ 𝔪^{i+1}` の**ノルム版**である。
⇒ ★**`h ∈ G_i ⟺ i ≤ t`**。前波の「1 つずれる」は解消した。 -/
theorem mem_ramification_iff [FiniteDimensional K M] {π : M} {n t i : ℕ}
    (h : M ≃ₐ[K] M) (hiso : ∀ w : M, ‖h w‖ = ‖w‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hn : Module.finrank K M = n)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    (hbr : ‖h π - π‖ = ‖π‖ ^ (t + 1)) :
    (∀ z : M, ‖z‖ ≤ 1 → ‖h z - z‖ ≤ ‖π‖ ^ (i + 1)) ↔ i ≤ t := by
  constructor
  · intro hall
    have hpi := hall π (le_of_lt hπ1)
    rw [hbr] at hpi
    by_contra hcon
    rw [not_le] at hcon
    have hlt : t + 1 < i + 1 := by omega
    have hstrict : ‖π‖ ^ (i + 1) < ‖π‖ ^ (t + 1) :=
      pow_lt_pow_right_of_lt_one₀ hπ0 hπ1 hlt
    linarith
  · intro hit z hz
    refine le_trans (norm_sub_apply_le_of_norm_le_one h hiso hπ0 hπ1 hn hvalK hbr hz) ?_
    exact pow_le_pow_of_le_one (le_of_lt hπ0) (le_of_lt hπ1) (by omega)

omit [IsUltrametricDist M] in
/-- ★★★**係数も整** —— `‖c·π^l‖ ≤ 1`(`l < n`)なら `‖c‖ ≤ 1`。

`‖c·π^l‖ = ‖π‖^{n·m + l} ≤ 1` から `n·m ≥ −l > −n`、ゆえに `m ≥ 0`、ゆえに `‖c‖ = ‖π‖^{n·m} ≤ 1`。
★これが `𝒪_M = 𝒪_K[π]`(単項生成)のノルム版の中身である。 -/
theorem norm_algebraMap_le_one_of_le_one {π : M} {n : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((n : ℤ) * m))
    {c : K} {l : ℕ} (hln : l < n) (hle : ‖algebraMap K M c * π ^ l‖ ≤ 1) :
    ‖algebraMap K M c‖ ≤ 1 := by
  rcases eq_or_ne c 0 with hc0 | hc0
  · rw [hc0, map_zero, norm_zero]; norm_num
  · obtain ⟨m, hm⟩ := hvalK c hc0
    have hval : ‖algebraMap K M c * π ^ l‖ = ‖π‖ ^ ((n : ℤ) * m + (l : ℤ)) := by
      rw [norm_mul, hm, norm_pow, ← zpow_natCast ‖π‖ l, ← zpow_add₀ (ne_of_gt hπ0)]
    have hE0 : 0 ≤ (n : ℤ) * m + (l : ℤ) := by
      by_contra hcon
      rw [not_le] at hcon
      have hgt : (1:ℝ) < ‖π‖ ^ ((n : ℤ) * m + (l : ℤ)) := one_lt_zpow_of_neg₀ hπ0 hπ1 hcon
      rw [← hval] at hgt
      linarith
    have hln' : (l : ℤ) < (n : ℤ) := by exact_mod_cast hln
    have hm0 : 0 ≤ m := by nlinarith [hE0, hln']
    rw [hm]
    calc ‖π‖ ^ ((n : ℤ) * m) ≤ ‖π‖ ^ (0 : ℤ) :=
          zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) (by positivity)
      _ = 1 := zpow_zero _


end Core
/-! ## `.src` と 公理 -/

def norm_sub_apply_le_of_norm_le_one.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def mem_ramification_iff.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

#print axioms norm_le_norm_pi_of_le_one
#print axioms norm_sub_apply_le_of_norm_le_one
#print axioms mem_ramification_iff
#print axioms norm_algebraMap_le_one_of_le_one

end RamNormBridge

end ABC3.Found.PGC
