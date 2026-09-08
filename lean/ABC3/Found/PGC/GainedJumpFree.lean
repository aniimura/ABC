import ABC3.Found.PGC.GainedJumpSeq

/-!
# [pGC] ★★★★★★`hu`(Hasse–Arf の中身)を落とす —— 0 の層は「縮小率 1」で通る

`Found/PGC/GainedJumpSeq.lean` の
`exists_norm_sub_algebraMap_le_prod_axDecay_of_conj` に残っていた仮説のうち、
★`hu : ∀ m ≤ k, p^m ≤ u m` **1 本だけ**が
「下付き跳びが `p` 冪で伸びる」= **Hasse–Arf の中身**であった。

★★本ファイルはその `hu` を落とす。★代わりに入るのは
`hi : 0 < i`(頂点が暴分岐)**1 本だけ**である。

★★★さらに §4 の最後で、★**跳びの列 `u` と `hjump` が丸ごと消える**
(`exists_norm_sub_algebraMap_le_prod_axDecay_of_conj_pure`)。
`u_m = if m = k then i else 0` を代入するだけで通るからで、これで仮説は

| 版 | 仮説 |
|---|---|
| `GainedJumpSeq...of_conj` | `he hu hjump hti hnormp heM hval hdeg htop hsep hsp hconj hne`(13) |
| 本ファイル `...of_conj_free` | `he hi hjump hti hnormp heM hval hdeg htop hsep hsp hconj hne`(13) |
| ★本ファイル `...of_conj_pure` | `he hi hnormp heM hval hdeg htop hsep hsp hconj hne`(★**11**) |

★★★**下付き分岐群の跳びの列は、この降下では 1 つも要らない。**
残る `i` は「全共役までの距離」(`hconj`)であって、跳びの**列**ではない。

## ★★★なぜ落ちるか —— 「1 未満に潰れた層は縮小率 1 で通す」

`GainedJumpSeq` の構成は `t_{j+1} = min(u_j, ⌊t_{j+2}/p⌋)`(上から `p` で割る鎖)で、
★`hu` はこの鎖が **`1` 以上に留まる**ためだけに使われていた。
鎖が `1` を切るのは `p^{k−j} > u_k` のとき、つまり
★**塔が深いのに頂点の跳びが小さい**ときである。

そこで鎖を `t_{j+1} := max(0, chain)` と切り上げる。すると

| 条件 | 0 に潰れた層で成り立つか |
|---|---|
| `ht0 : 0 ≤ t_{j+1}` | ★定義から |
| `hstepZ : p·t_{j+1} ≤ t_{j+2}` | ★`p·0 = 0 ≤ max(0,·)` |
| `hlayerZ : (p−1)t_{j+1} ≤ p^{j+1}e` | ★`0 ≤ p^{j+1}e` |
| `hti : i = t_{k+1}` | ★頂点は `max(0,u_k) = u_k`(`hti` と `i : ℕ` から `u_k ≥ 0`) |
| ★`hstep : ‖τ^{p^j}z − z‖ ≤ ‖π‖^{t_{j+1}}‖z‖` | ★★**これだけが残る問題**(`‖π‖^0 = 1`) |

★★最後の 1 行が本ファイルの数学である。`t_{j+1} = 0` の層では
`‖τ^{p^j}z − z‖ ≤ ‖z‖`(縮小率 1)を示せばよく、これは

```
‖σπ‖ ≤ max(‖σπ − π‖, ‖π‖) = ‖π‖   ⟹   ‖σz‖ ≤ ‖z‖   ⟹   ‖σz − z‖ ≤ max(‖σz‖,‖z‖) = ‖z‖
```

で出る。真ん中は**桁の分離**(`CyclicJumpNorm.nnnorm_sum_digit_eq_sup`、木に在った)から:
`σ(Σ a_jπ^j) = Σ a_j(σπ)^j` の各桁が `‖a_jπ^j‖` 以下、和は `sup` 以下、`sup = ‖z‖`。
★必要な入力は `‖σπ − π‖ ≤ ‖π‖`(`hjump0`)**だけ**である。

## ★★★`hjump0` は `hconj` から**無料で**出る(★これで新しい仮説が 0 になる)

`hjump0 : ∀ j < k, ‖τ^{p^j}π − π‖ ≤ ‖π‖` は
`exists_norm_sub_algebraMap_le_prod_axDecay_of_jumps_free` では仮説だが
(`hjump` と `0 ≤ u_j` から出る)、
★`exists_norm_sub_algebraMap_le_prod_axDecay_of_conj_free` では**導出できる**:

`τ^{p^j}π` は `minpoly F π` の根なので、`hconj` により
`τ^{p^j}π ≠ π` なら `‖τ^{p^j}π − π‖ = ‖π‖^{i+1} ≤ ‖π‖`(`i ≥ 0`、`‖π‖ < 1`)、
`τ^{p^j}π = π` なら `0 ≤ ‖π‖`。★どちらの枝も `hconj` だけで閉じる。

★★`GainedBridgeSupply.norm_algHom_sub_le_mul_of_break` は縮小率に `θ < 1` を要求する
(`‖σπ−π‖ < ‖π‖` から桁の**等式** `norm_digitSum_sub_mul_norm_eq` を使うため)。
★θ = 1 の端では等式が壊れるので、本ファイルは**不等式だけの別証**を置いた
(`norm_algHom_sub_le_self`)。★これが「落ちなかった 1 本」を落とした鍵である。

## ★★★予算は変わらない(★数値で確かめた)

`t` が 0 に潰れても出口の定数 `∏_{j∈[1,k+1]} axDecay p j` は**同じ**である。
`GainedTowerDescent.gainedLoss` は `max` で定義されているので `Λ ≥ 0` は保たれ、
`gainedLoss_fits`(`(p−1)²Λ ≤ (p^{k+1}−p)e`)は `ht0`/`hstepZ`/`hlayerZ` しか使わない。

★退化した例で勘定を閉じた(`Numeric.free_zero_low`、機械検算):
`p = 3`, `k = 1`, `e = 1`, 下付き跳び `u = (u₀,u₁) = (0,4)`。
本ファイルの構成は `t = (0, 4, 12, …)` を作り、

* `hstepZ`: `3·0 = 0 ≤ 4` ✓、`3·4 = 12 ≤ 12` ✓
* `hlayerZ`: `2·0 = 0 ≤ 3·1` ✓、`2·4 = 8 ≤ 9·1` ✓
* `jumpSum t = (0, 0, 4)`、`Λ₂ = max(4 − 2·0, max(0, 4−0−0) + 3·0) = 4`
* `gainedLoss_fits`(`k = 2`): `(3−1)²·4 = 16 ≤ (3³−3)·1 = 24` ✓

★`GainedJumpSeq.exists_jump_seq` の `hu` は `3⁰ = 1 ≤ u₀ = 0` を要求するので
**この例は扱えなかった**。★予算(`gainedLoss_fits`)には余裕があり、
落ちていたのは数学ではなく構成であった。

## ★★★測定(コマンドつき)

| 部品 | 在庫 | 測定 |
|---|---|---|
| 桁の分離 `‖Σa_jπ^j‖ = sup‖a_jπ^j‖` | ★★**木に在った** | `grep -n "nnnorm_sum_digit_eq_sup" lean/ABC3/Found/PGC/CyclicJumpNorm.lean`(310 行) |
| `σ(Σa_jπ^j) = Σa_j(σπ)^j` | ★木に在った | 同 301 行 `map_digitSum`(★`σ` に等長性も全単射性も要らない) |
| 超距離の和の上界 | mathlib に在った | `#check @IsUltrametricDist.nnnorm_sum_le_of_forall_le` |
| `‖x − y‖ ≤ max‖x‖‖y‖` | ★**名前は無い** | `IsUltrametricDist.norm_sub_le_max` は `Unknown constant`。`norm_add_le_max` に `sub_eq_add_neg` で落とす |
| Hasse–Arf の橋 | ★**要らなくなった** | 本ファイルは `HasseArf*.lean` を import しない |

## ★★★閉じていないもの(★正確に)

1. ★`hi : 0 < i` は落ちない。`i = 0` は「頂点が順分岐」で、そこでは
   `‖τ^{p^k}π − π‖ = ‖π‖` となり、降下の要である桁の等式
   (`CyclicJumpNorm.norm_digitSum_sub_mul_norm_eq` の仮定 `‖u − π‖ < ‖π‖`)が**壊れる**。
   ★これは構成の問題ではなく境界そのもので、`hbreak` / `hti` と同じ設定の一部である。
2. ★`hjump` は `...of_conj_pure` で**落ちた**。★`GainedJumpSeq` の記述
   「`hjump` は下付き分岐群の定義そのもので消えない」は**外れ**である
   (下の層の跳びは結論の定数に効かない ——
   `gainedLoss_fits` の予算が最悪の場合でも足りる)。
3. (B) の配管 `hnormp` / `heM` / `hval` / `hdeg` / `htop` / `hsep` / `hsp` /
   `hconj` / `hne` は本ファイルでも触っていない。★これは新しい `structure` を作る仕事で、
   `Skeleton/PGC/Setup.lean:40` の `PAdicLocalField` は
   `carrier / Field / Algebra ℚ_p / FiniteDimensional` の 4 つしか持たない
   (★素元もノルムも分岐指数も無いので、そこからは `hval` も `heM` も出ない)。

## 逸脱の記録(CLAUDE.md「逸脱」)

1. `GainedJumpSeq` / `GainedBridgeSupply` / `CyclicJumpNorm` は**読むだけ**で 1 行も
   書き換えていない。弱めた版はすべて本ファイルに別名で置いた。
2. ★`GainedJumpSeq.exists_jump_seq` は `1 ≤ t_{j+1}` を出すが、本ファイルの
   `exists_jump_seq_free` は `0 ≤ t_{j+1}` と「`t_{j+1} = 0` か `t_{j+1} ≤ u_j`」の
   **選言**を出す。★選言にしたのは、0 に潰れた層で `u_j` を経由しないためである。
3. 原典(Serre *Corps Locaux* IV)は付値言語。本ファイルはノルム言語
   (`v_M(z) = m` は `‖z‖ = ‖π‖^m`)。`GainedJumpSeq` と同じ約束である。
-/
open Finset

namespace ABC3.Found.PGC

namespace GainedJumpFree

/-! ## §1 抽象核(`ℤ` のみ) —— 0 で切り上げた鎖 -/

/-- ★★★**抽象核**: `GainedJumpSeq.exists_jump_seq` から
`hu : ∀ m ≤ k, p^m ≤ u m` を落とした形。

★入力は `0 ≤ e`、`0 ≤ u k`、頂点の上界 `(p−1)u_k ≤ p^{k+1}e` の 3 本だけである。
★出力の `t` は `1 ≤ t_{j+1}` を落とし、代わりに
「`t_{j+1} = 0` か `t_{j+1} ≤ u_j`」という**選言**を出す。

★分岐・付値・Galois の語彙が 1 語も出ない。 -/
theorem exists_jump_seq_free {p : ℤ} (hp : 2 ≤ p) {u : ℕ → ℤ} {e : ℤ} (k : ℕ)
    (he : 0 ≤ e) (huk : 0 ≤ u k)
    (hbnd : (p - 1) * u k ≤ p ^ (k + 1) * e) :
    ∃ t : ℕ → ℤ,
      t (k + 1) = u k ∧
      (∀ j, 0 ≤ t (j + 1)) ∧
      (∀ j, j ≤ k → t (j + 1) = 0 ∨ t (j + 1) ≤ u j) ∧
      (∀ j, p * t (j + 1) ≤ t (j + 2)) ∧
      (∀ j, (p - 1) * t (j + 1) ≤ p ^ (j + 1) * e) := by
  have hp0 : (0 : ℤ) < p := by linarith
  refine ⟨fun n => max 0 (GainedJumpSeq.seq p u k (n - 1)), ?_, ?_, ?_, ?_, ?_⟩
  · show max 0 (GainedJumpSeq.seq p u k (k + 1 - 1)) = u k
    have he1 : k + 1 - 1 = k := by omega
    rw [he1, GainedJumpSeq.seq_top p u k]
    exact max_eq_right huk
  · intro j; simp only [Nat.add_sub_cancel]; exact le_max_left _ _
  · intro j hj
    simp only [Nat.add_sub_cancel]
    rcases le_or_gt (GainedJumpSeq.seq p u k j) 0 with h | h
    · exact Or.inl (max_eq_left h)
    · exact Or.inr (by rw [max_eq_right (le_of_lt h)]; exact GainedJumpSeq.seq_le_of_le p u hj)
  · intro j
    have e1 : j + 1 - 1 = j := by omega
    have e2 : j + 2 - 1 = j + 1 := by omega
    simp only [e1, e2]
    have hstep := GainedJumpSeq.mul_seq_le_seq_succ hp u k j
    rcases le_or_gt (GainedJumpSeq.seq p u k j) 0 with h | h
    · rw [max_eq_left h]
      simp
    · rw [max_eq_right (le_of_lt h)]
      exact le_trans hstep (le_max_right _ _)
  · intro j
    simp only [Nat.add_sub_cancel]
    rcases le_or_gt (GainedJumpSeq.seq p u k j) 0 with h | h
    · rw [max_eq_left h]
      have : (0 : ℤ) ≤ p ^ (j + 1) * e := mul_nonneg (le_of_lt (pow_pos hp0 _)) he
      linarith
    · rw [max_eq_right (le_of_lt h)]
      exact GainedJumpSeq.sub_one_mul_seq_le hp hbnd j

/-! ## §2 抽象核(ノルム) —— ★縮小率 `1` の端 -/

section NormCore

variable {F M : Type*} [Field F] [NormedField M] [IsUltrametricDist M] [Algebra F M]

omit [Algebra F M] in
/-- 超距離の差の上界。★mathlib に `IsUltrametricDist.norm_sub_le_max` は**無い**
(`Unknown constant`)ので `norm_add_le_max` から作る。 -/
theorem norm_sub_le_max' (x y : M) : ‖x - y‖ ≤ max ‖x‖ ‖y‖ := by
  simpa [sub_eq_add_neg] using IsUltrametricDist.norm_add_le_max x (-y)

/-- ★★★**縮小率 `1` の端**: `‖σπ − π‖ ≤ ‖π‖` だけから `‖σ z − z‖ ≤ ‖z‖`。

★`GainedBridgeSupply.norm_algHom_sub_le_mul_of_break` は `θ < 1` を要求する
(桁の**等式** `norm_digitSum_sub_mul_norm_eq` が `‖σπ−π‖ < ‖π‖` を要るため)。
本補題は等式を使わず、**桁の分離の `sup` と超距離の和の上界だけ**で通す。

証明は 3 行:
`‖σπ‖ ≤ max(‖σπ−π‖,‖π‖) = ‖π‖` → 各桁 `‖a_j(σπ)^j‖ ≤ ‖a_jπ^j‖ ≤ ‖z‖` →
`‖σz‖ ≤ ‖z‖` → `‖σz − z‖ ≤ max(‖σz‖,‖z‖) = ‖z‖`。 -/
theorem norm_algHom_sub_le_self {p : ℕ} (hp : 0 < p) {π : M} (σ : M →ₐ[F] M)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1)
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hbrk : ‖σ π - π‖ ≤ ‖π‖)
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) (z : M) :
    ‖σ z - z‖ ≤ ‖z‖ := by
  have hint : IsIntegral F π := by
    by_contra hc
    rw [minpoly.eq_zero hc, Polynomial.natDegree_zero] at hdeg
    omega
  obtain ⟨a, rfl⟩ := exists_digitSum_of_adjoin_eq_top hint hdeg htop z
  have hσπ : ‖σ π‖ ≤ ‖π‖ := by
    have h := norm_sub_le_max' (σ π - π) (-π)
    rw [sub_neg_eq_add, sub_add_cancel, norm_neg] at h
    exact le_trans h (max_le hbrk le_rfl)
  have hσπnn : ‖σ π‖₊ ≤ ‖π‖₊ := by exact_mod_cast hσπ
  have hsup : ‖digitSum p π a‖₊
      = (Finset.range p).sup fun j => ‖algebraMap F M (a j) * π ^ j‖₊ := by
    rw [digitSum]
    exact nnnorm_sum_digit_eq_sup hπ0 hπ1 hval a fun j hj => Finset.mem_range.mp hj
  have hσle : ‖σ (digitSum p π a)‖₊ ≤ ‖digitSum p π a‖₊ := by
    rw [map_digitSum, digitSum]
    refine IsUltrametricDist.nnnorm_sum_le_of_forall_le fun j hj => ?_
    rw [hsup]
    refine le_trans ?_ (Finset.le_sup hj)
    simp only [nnnorm_mul, nnnorm_pow]
    gcongr
  have hσle' : ‖σ (digitSum p π a)‖ ≤ ‖digitSum p π a‖ := by exact_mod_cast hσle
  exact le_trans (norm_sub_le_max' _ _) (max_le hσle' le_rfl)

/-- ★★`t_{j+1} ≥ 0` の**すべての**層で使える形にまとめたもの。

`c = 0` の層は `norm_algHom_sub_le_self`(縮小率 1)、
`c ≥ 1` の層は `GainedBridgeSupply.norm_algHom_sub_le_mul_of_break`(縮小率 `‖π‖^c < 1`)。
★これが `GainedJumpSeq` の `ht1 : 1 ≤ t_{j+1}` を `ht0 : 0 ≤ t_{j+1}` に落とす橋である。 -/
theorem norm_algHom_sub_le_zpow_mul {p : ℕ} (hp : 0 < p) {π : M} (σ : M →ₐ[F] M) {c : ℤ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hc0 : 0 ≤ c)
    (hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1)
    (hval : ∀ d : F, d ≠ 0 → ∃ m : ℤ, ‖algebraMap F M d‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hbrk : ‖σ π - π‖ ≤ ‖π‖ ^ (c + 1))
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) (z : M) :
    ‖σ z - z‖ ≤ ‖π‖ ^ c * ‖z‖ := by
  rcases eq_or_lt_of_le hc0 with hc | hc
  · rw [← hc, zpow_zero, one_mul]
    refine norm_algHom_sub_le_self hp σ hπ0 hπ1 hval ?_ hdeg htop z
    rw [← hc] at hbrk
    simpa using hbrk
  · refine GainedBridgeSupply.norm_algHom_sub_le_mul_of_break hp σ hπ0 hπ1
      (zpow_lt_one₀ hπ0 hπ1 hc) hchar hval ?_ hdeg htop z
    rwa [zpow_add_one₀ (ne_of_gt hπ0)] at hbrk

end NormCore

/-! ## §3 出口 —— ★`hu`(Hasse–Arf)が `hjump0`(`G_0`)と `hi`(暴分岐)に落ちる -/

section Exit

variable {F M : Type*} [Field F] [NormedField M] [IsUltrametricDist M] [Algebra F M]

/-- ★★★★★★**`GainedJumpSeq.exists_norm_sub_algebraMap_le_prod_axDecay_of_jumps` から
`hu : ∀ m ≤ k, p^m ≤ u m` を落とした形**。

代わりに入るのは
`hjump0 : ∀ j < k, ‖τ^{p^j}π − π‖ ≤ ‖π‖`(★`τ^{p^j} ∈ G_0`。`hjump` と `0 ≤ u_j` から出る)と
`hi : 0 < i`(★頂点が暴分岐)の 2 本で、★どちらも Hasse–Arf を使わない。
★§4 ではこの `hjump0` が `hconj` から**無料で**出るので、増えるのは `hi` の 1 本だけになる。

★★`hu` は「跳びが `p` 冪で伸びる」= Hasse–Arf の中身だった。
本形は跳びの列が `1` を切る層を `t_{j+1} = 0`(縮小率 1)で通すので、
★**跳びの伸びを一切仮定しない**。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_jumps_free
    {p : ℕ} [Fact p.Prime] {e i k : ℕ} {π : M} {u : ℕ → ℤ}
    (τ : M →ₐ[F] M) (he : 0 < e) (hi : 0 < i)
    (hjump0 : ∀ j, j < k → ‖(τ ^ p ^ j) π - π‖ ≤ ‖π‖)
    (hjump : ∀ j, j < k → ‖(τ ^ p ^ j) π - π‖ ≤ ‖π‖ ^ (u j + 1))
    (hbnd : ((p : ℤ) - 1) * u k ≤ (p : ℤ) ^ (k + 1) * (e : ℤ))
    (hti : (i : ℤ) = u k)
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (hbreak : ‖(τ ^ p ^ k) π - π‖ = ‖π‖ ^ (i + 1))
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤) (x : M) :
    ∃ y : F, ‖x - algebraMap F M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ := by
  have hp' : p.Prime := Fact.out
  have hpZ : (2 : ℤ) ≤ (p : ℤ) := by exact_mod_cast hp'.two_le
  have huk : (0 : ℤ) ≤ u k := by rw [← hti]; exact Int.natCast_nonneg i
  obtain ⟨t, htk, ht0, htu, hstepZ, hlayerZ⟩ :=
    exists_jump_seq_free hpZ (e := (e : ℤ)) k (Int.natCast_nonneg e) huk hbnd
  have hti' : (i : ℤ) = t (k + 1) := by rw [htk]; exact hti
  have hpR : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp'.one_lt
  have hpR0 : (0 : ℝ) < (p : ℝ) := lt_trans zero_lt_one hpR
  have hEne : p ^ (k + 1) * e ≠ 0 := Nat.mul_ne_zero (pow_ne_zero _ hp'.ne_zero) he.ne'
  have hπpow : ‖π‖ ^ (p ^ (k + 1) * e) = ((p : ℝ))⁻¹ := by rw [← heM]; exact hnormp
  have hplt : ‖(p : M)‖ < 1 := by rw [hnormp, inv_lt_one_iff₀]; exact Or.inr hpR
  have hπ0 : 0 < ‖π‖ := by
    rcases lt_or_eq_of_le (norm_nonneg π) with h | h
    · exact h
    · exfalso
      rw [← h, zero_pow hEne] at hπpow
      exact absurd hπpow.symm (ne_of_gt (inv_pos.mpr hpR0))
  have hπ1 : ‖π‖ < 1 := by
    refine GainedBridgeSupply.lt_one_of_pow_lt_one (E := p ^ (k + 1) * e) ?_
    rw [hπpow, inv_lt_one_iff₀]; exact Or.inr hpR
  have hchar : ∀ j, 0 < j → j < p → ‖(j : M)‖ = 1 :=
    fun j hj0 hjp => norm_natCast_eq_one_of_lt_prime hp' hplt hj0 hjp
  set τ' : M →+ M := AddMonoidHom.mk' (⇑τ) (fun a b => map_add τ a b) with hτ'
  have hcoe : ⇑τ' = ⇑τ := rfl
  have hpowcoe : ∀ j : ℕ, (⇑τ')^[p ^ j] = ⇑(τ ^ p ^ j) := by
    intro j; rw [hcoe, AlgHom.coe_pow]
  have hstep : ∀ j, j < k → ∀ z : M, ‖(⇑τ')^[p ^ j] z - z‖ ≤ ‖π‖ ^ (t (j + 1)) * ‖z‖ := by
    intro j hj z
    rw [hpowcoe j]
    refine norm_algHom_sub_le_zpow_mul hp'.pos (τ ^ p ^ j) hπ0 hπ1 (ht0 j) hchar hval ?_
      hdeg htop z
    rcases htu j (le_of_lt hj) with h | h
    · rw [h]
      simpa using hjump0 j hj
    · refine le_trans (hjump j hj) ?_
      refine zpow_le_zpow_right_of_le_one₀ hπ0 (le_of_lt hπ1) ?_
      omega
  have hσ : ∀ z : M, (τ ^ p ^ k) z = (⇑τ')^[p ^ k] z := by intro z; rw [hpowcoe k]
  obtain ⟨y, hy⟩ := GainedJumpSeq.exists_norm_sub_algebraMap_le_prod_axDecay_lt
    (F := F) (M := M) (p := p) (e := e) (i := i) (k := k) (π := π) τ' (τ ^ p ^ k) t he hσ ht0
    hstep hstepZ hlayerZ hi hti' hnormp heM hbreak hval hdeg htop x
  exact ⟨y, hy⟩

end Exit

/-! ## §4 出口 —— `hbnd` と `hbreak` も落とした形(`GainedJumpSeq` §9 の `hu`-free 版) -/

section Exit2

variable {F M : Type*} [Field F] [NormedField M] [IsUltrametricDist M] [Algebra F M]

omit [IsUltrametricDist M] in
/-- ★★`hconj`(全共役の距離)から `hjump0`(縮小率 1 の入力)が**無料で**出る。

`σ π` は `minpoly F π` の根なので、動くなら差はちょうど `‖π‖^{i+1}`、
`‖π‖ < 1` かつ `i ≥ 0` なので `≤ ‖π‖`。動かないなら `0 ≤ ‖π‖`。 -/
theorem norm_algHom_sub_le_of_conj {i : ℕ} {π : M} (hπ1 : ‖π‖ < 1) (σ : M →ₐ[F] M)
    (hconj : ∀ a : M, ((minpoly F π).map (algebraMap F M)).IsRoot a → a ≠ π →
      ‖π - a‖ = ‖π‖ ^ (i + 1)) : ‖σ π - π‖ ≤ ‖π‖ := by
  rcases eq_or_ne (σ π) π with h | h
  · simp [h, norm_nonneg]
  · rw [GainedJumpSeq.norm_algHom_sub_eq_of_conj σ hconj h]
    calc ‖π‖ ^ (i + 1) ≤ ‖π‖ ^ 1 :=
          pow_le_pow_of_le_one (norm_nonneg π) (le_of_lt hπ1) (by omega)
      _ = ‖π‖ := pow_one _

/-- ★★★★★★**`AxLemmaGraded` の形。Hasse–Arf を使わない版**。

`GainedJumpSeq.exists_norm_sub_algebraMap_le_prod_axDecay_of_conj` と同じ結論を、
★`hu : ∀ m ≤ k, p^m ≤ u m`(**Hasse–Arf の中身**)**なし**で出す。

残る仮説は 13 本:
`he`(`0 < e`)、★`hi`(`0 < i`、頂点が暴分岐)、
`hjump`(下付き分岐群の定義)、`hti`(`i` の定義)、
`hnormp` / `heM`(正規化)、`hval`(値群)、
`hdeg` / `htop` / `hsep` / `hsp` / `hconj` / `hne`(`M = F(π)` が `p` 次 Galois)。

★★`GainedJumpSeq` 版と**本数は同じ**で、`hu`(Hasse–Arf)が `hi` に置き換わっている。
★★**跳びの列の伸びに関する仮定は 1 つも残っていない。** -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_conj_free
    {p : ℕ} [Fact p.Prime] {e i k : ℕ} {π : M} {u : ℕ → ℤ}
    (τ : M →ₐ[F] M) (he : 0 < e) (hi : 0 < i)
    (hjump : ∀ j, j < k → ‖(τ ^ p ^ j) π - π‖ ≤ ‖π‖ ^ (u j + 1))
    (hti : (i : ℤ) = u k)
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤)
    (hsep : (minpoly F π).Separable)
    (hsp : ((minpoly F π).map (algebraMap F M)).Splits)
    (hconj : ∀ a : M, ((minpoly F π).map (algebraMap F M)).IsRoot a → a ≠ π →
      ‖π - a‖ = ‖π‖ ^ (i + 1))
    (hne : (τ ^ p ^ k) π ≠ π) (x : M) :
    ∃ y : F, ‖x - algebraMap F M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ := by
  have hp' : p.Prime := Fact.out
  have hp1 : 1 ≤ p := hp'.one_lt.le
  have hpR : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp'.one_lt
  have hpR0 : (0 : ℝ) < (p : ℝ) := lt_trans zero_lt_one hpR
  have hEne : p ^ (k + 1) * e ≠ 0 := Nat.mul_ne_zero (pow_ne_zero _ hp'.ne_zero) he.ne'
  have hπpow : ‖π‖ ^ (p ^ (k + 1) * e) = ((p : ℝ))⁻¹ := by rw [← heM]; exact hnormp
  have hπ0 : 0 < ‖π‖ := by
    rcases lt_or_eq_of_le (norm_nonneg π) with h | h
    · exact h
    · exfalso
      rw [← h, zero_pow hEne] at hπpow
      exact absurd hπpow.symm (ne_of_gt (inv_pos.mpr hpR0))
  have hπ1 : ‖π‖ < 1 := by
    refine GainedBridgeSupply.lt_one_of_pow_lt_one (E := p ^ (k + 1) * e) ?_
    rw [hπpow, inv_lt_one_iff₀]; exact Or.inr hpR
  have hint : IsIntegral F π := by
    by_contra hc
    rw [minpoly.eq_zero hc, Polynomial.natDegree_zero] at hdeg
    exact hp'.ne_zero hdeg.symm
  have hbreak : ‖(τ ^ p ^ k) π - π‖ = ‖π‖ ^ (i + 1) :=
    GainedJumpSeq.norm_algHom_sub_eq_of_conj (τ ^ p ^ k) hconj hne
  -- ★`hjump0`(`t_{j+1} = 0` の層で使う縮小率 1)は `hconj` から**無料で**出る:
  -- `τ^{p^j}π` は最小多項式の根なので、動くなら差はちょうど `‖π‖^{i+1} ≤ ‖π‖`。
  have hjump0 : ∀ j, j < k → ‖(τ ^ p ^ j) π - π‖ ≤ ‖π‖ :=
    fun j _ => norm_algHom_sub_le_of_conj hπ1 (τ ^ p ^ j) hconj
  have hnat : (p - 1) * i ≤ p ^ (k + 1) * e :=
    GainedJumpSeq.sub_one_mul_jump_le_of_splits hp'.pos hπ0 hπ1 hval hint hdeg hsep hsp hconj heM
  have hcast : ((p - 1 : ℕ) : ℤ) = (p : ℤ) - 1 := by
    push_cast [Nat.cast_sub hp1]; ring
  have hbnd : ((p : ℤ) - 1) * u k ≤ (p : ℤ) ^ (k + 1) * (e : ℤ) := by
    have h3 : (((p - 1) * i : ℕ) : ℤ) ≤ ((p ^ (k + 1) * e : ℕ) : ℤ) := by exact_mod_cast hnat
    rw [Nat.cast_mul, hcast, Nat.cast_mul, Nat.cast_pow] at h3
    rw [← hti]
    exact h3
  exact exists_norm_sub_algebraMap_le_prod_axDecay_of_jumps_free τ he hi hjump0 hjump hbnd hti
    hnormp heM hbreak hval hdeg htop x

/-- ★★★★★★**跳びの列 `u` そのものが消えた形**(★本ファイルの最終出口)。

`exists_norm_sub_algebraMap_le_prod_axDecay_of_conj_free` の `u` に
★`u_m = if m = k then i else 0`(**下の層は全部 0**)を代入しただけである。

★★これが通るのは、結論の定数 `∏_{j∈[1,k+1]} axDecay p j` が
「下の層で 1 つも得をしない」最悪の場合でも足りるからで、
`gainedLoss_fits`(`(p−1)²Λ ≤ (p^{k+2}−p)e`)が `Λ = i` と
`hlayerZ`(`(p−1)i ≤ p^{k+1}e`)から
`(p−1)²i ≤ (p−1)p^{k+1}e = (p^{k+2}−p^{k+1})e ≤ (p^{k+2}−p)e` で閉じる。

★★★**下付き分岐群の跳びの列は、この降下では 1 つも要らない。**
残るのは `i`(頂点 = 全共役の距離、`hconj` が決める)だけである。

残る仮説は 11 本:
`he`(`0 < e`)、`hi`(`0 < i`、頂点が暴分岐)、
`hnormp` / `heM`(正規化)、`hval`(値群)、
`hdeg` / `htop` / `hsep` / `hsp` / `hconj` / `hne`(`M = F(π)` が `p` 次 Galois)。 -/
theorem exists_norm_sub_algebraMap_le_prod_axDecay_of_conj_pure
    {p : ℕ} [Fact p.Prime] {e i k : ℕ} {π : M}
    (τ : M →ₐ[F] M) (he : 0 < e) (hi : 0 < i)
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹)
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p ^ (k + 1) * e))
    (hval : ∀ c : F, c ≠ 0 → ∃ m : ℤ, ‖algebraMap F M c‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hdeg : (minpoly F π).natDegree = p)
    (htop : Algebra.adjoin F ({π} : Set M) = ⊤)
    (hsep : (minpoly F π).Separable)
    (hsp : ((minpoly F π).map (algebraMap F M)).Splits)
    (hconj : ∀ a : M, ((minpoly F π).map (algebraMap F M)).IsRoot a → a ≠ π →
      ‖π - a‖ = ‖π‖ ^ (i + 1))
    (hne : (τ ^ p ^ k) π ≠ π) (x : M) :
    ∃ y : F, ‖x - algebraMap F M y‖
      ≤ (∏ j ∈ Finset.Icc 1 (k + 1), axDecay p j) * ‖τ x - x‖ := by
  have hp' : p.Prime := Fact.out
  have hpR : (1 : ℝ) < (p : ℝ) := by exact_mod_cast hp'.one_lt
  have hpR0 : (0 : ℝ) < (p : ℝ) := lt_trans zero_lt_one hpR
  have hEne : p ^ (k + 1) * e ≠ 0 := Nat.mul_ne_zero (pow_ne_zero _ hp'.ne_zero) he.ne'
  have hπpow : ‖π‖ ^ (p ^ (k + 1) * e) = ((p : ℝ))⁻¹ := by rw [← heM]; exact hnormp
  have hπ1 : ‖π‖ < 1 := by
    refine GainedBridgeSupply.lt_one_of_pow_lt_one (E := p ^ (k + 1) * e) ?_
    rw [hπpow, inv_lt_one_iff₀]; exact Or.inr hpR
  refine exists_norm_sub_algebraMap_le_prod_axDecay_of_conj_free
    (u := fun m => if m = k then (i : ℤ) else 0) τ he hi ?_ ?_ hnormp heM hval hdeg htop
    hsep hsp hconj hne x
  · intro j hj
    have hzero : (if j = k then (i : ℤ) else 0) = 0 := if_neg (by omega)
    rw [hzero]
    simpa using norm_algHom_sub_le_of_conj hπ1 (τ ^ p ^ j) hconj
  · simp

end Exit2

/-! ## §5 数値の検算 —— ★`GainedJumpSeq` が扱えなかった退化例 -/

namespace Numeric

/-- ★`p = 3`, `k = 1`, 下付き跳び `u = (u₀,u₁) = (0,4)`(★塔が深いのに下の跳びが 0)。

`GainedJumpSeq.exists_jump_seq` の `hu` は `3^0 = 1 ≤ u₀ = 0` を要求するので**通らない**。
★本ファイルの構成は `t = (max(0, min(0,⌊4/3⌋)), 4) = (0,4)` を作り、
`hstepZ`(`3·0 = 0 ≤ 4`)を満たす。`t₁ = 0` の層は縮小率 1 で通す。 -/
theorem free_zero_low :
    ¬ ((3 : ℤ) ^ 0 ≤ (fun m => if m = 0 then (0 : ℤ) else 4) 0) ∧
    max 0 (GainedJumpSeq.seq 3 (fun m => if m = 0 then (0 : ℤ) else 4) 1 0) = 0 ∧
    max 0 (GainedJumpSeq.seq 3 (fun m => if m = 0 then (0 : ℤ) else 4) 1 1) = 4 ∧
    3 * max 0 (GainedJumpSeq.seq 3 (fun m => if m = 0 then (0 : ℤ) else 4) 1 0)
      ≤ max 0 (GainedJumpSeq.seq 3 (fun m => if m = 0 then (0 : ℤ) else 4) 1 1) := by
  refine ⟨by norm_num, ?_, ?_, ?_⟩ <;>
    norm_num [GainedJumpSeq.seq, GainedJumpSeq.chain]

end Numeric

/-! ## §6 `.src`(原典の対応箇所) -/

def exists_jump_seq_free.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def norm_algHom_sub_le_self.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_prod_axDecay_of_jumps_free.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_prod_axDecay_of_conj_free.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_prod_axDecay_of_conj_pure.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §7 使っている公理の一覧 -/

#print axioms exists_jump_seq_free
#print axioms norm_sub_le_max'
#print axioms norm_algHom_sub_le_self
#print axioms norm_algHom_sub_le_zpow_mul
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_jumps_free
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_conj_free
#print axioms norm_algHom_sub_le_of_conj
#print axioms exists_norm_sub_algebraMap_le_prod_axDecay_of_conj_pure
#print axioms Numeric.free_zero_low

end GainedJumpFree

end ABC3.Found.PGC
