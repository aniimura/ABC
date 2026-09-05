import ABC3.Found.PGC.UnramifiedFrobenius
import ABC3.Found.Teichmuller

/-!
# `K^ur` のノルム 1 の元は `p ∤ n` 乗根を `K^ur` の中に持つ

経路 C(`ResearchPaper/pgc-goal.md`)のノード 1
`units_isNthPow_unramifiedClosure` ——

> `p ∤ n` なら `𝒪_{K^ur}^×` の元は `K^ur` の中で `n` 乗数である

の形式化。Kummer 理論を `K^ur` 上で回すとき、`𝒪_{K^ur}^×` が `n` 可除で
あることが `N_n(Γ_K)` の上界を出す入口になる。

## ★整数環を名指ししない書き方(設計上の要点)

ABC3 には **`𝒪_{K^ur}`(`unramifiedClosure K` の整数環)が存在しない**。
`unramifiedClosure K : IntermediateField K.carrier K.closure` はあるが、
整数環は**有限部分拡大ごとの `adjoinIntegers K z : Subring ↥K⟮z⟯`** しか無く、
`HenselianLocalRing` インスタンスもそちらにしか付いていない。

そこで **整数環を一切名指しせず `‖·‖ = 1` で書く**:

```
x ∈ unramifiedClosure K → ‖x‖ = 1 → ∃ y ∈ unramifiedClosure K, y ^ n = x
```

`‖·‖` は `Found/PGC/LocalFieldNorm.lean::closureNormedField`(スペクトル
ノルム)。`K^ur` は `K.closure` の中間体なので、`𝒪_{K^ur}^×` は
「`K^ur` の元で `‖·‖ = 1` のもの」と同じであり、上の形は主張として同値。
新しい部分環を構成しないので `adjoinIntegers` と `𝒪[(adjoinField K x).carrier]`
の境界(`tools/lean-idioms.md` #69)を一度も越えない。

## 証明の筋(★採った道)

`x ∈ K^ur` は `mem_unramifiedClosure_iff_isUnramified` により
**`K(x)/K` 自身が不分岐**ということなので、`z` を別に取る必要が無く、
`A := adjoinIntegers K x`(`K(x)` のノルム `≤ 1` の元)で全部が済む。

1. `‖x‖ = 1` なので `x` は `A` の単数(`isUnit_of_norm_eq_one`)。
2. `A` の剰余体は有限体 `𝔽_Q`(`Q = q^f`、したがって `p ∣ Q`、
   ゆえに **`p ∤ Q − 1`**)。
3. **Teichmüller 分解**(`Found/Teichmuller.lean::exists_teichmullerRep`、
   Hensel):`x = ζ · u`、`ζ^{Q−1} = 1`、`u` の剰余は `1`。
4. **主単数 `u` は `A` の中で `n` 乗根を持つ**(Hensel)——`X^n − u` の
   剰余は `X^n − 1` で、`1` が根、微分は `n ≠ 0`(`p ∤ n`)。
   これが `exists_pow_eq_of_residue_eq_one`。
5. **1 の冪根 `ζ` は `K^ur` の中で `n` 乗根を持つ**——`p ∤ n(Q−1)` なので
   `μ_{n(Q−1)} ⊆ K^ur`(在庫、第 1005 `exists_isPrimitiveRoot_mem_
   unramifiedClosure`)。原始 `n(Q−1)` 乗根 `ξ` を取れば `ξ^n` は原始
   `Q−1` 乗根なので `ζ = (ξ^n)^j = (ξ^j)^n`。
6. `y := ξ^j · v`。

★`pgc-goal.md` の当初案は「より大きな有限不分岐拡大へ上げ、剰余体
`𝔽_{q^M}` で `n` 乗根を取る」だったが、**塔を登る必要は無かった**:
`μ_{n(Q−1)} ⊆ K^ur` が既に在庫にあるので、`n` 乗根が要るのは
「1 の冪根の部分」だけで、そこは `K^ur` の中で直接取れる。主単数の
部分は `A` を出ずに Hensel で済む。`adjoin_le_of_dvd`(次数の整除で
塔に載せる)も `exists_isUnramifiedAdjoin K M` も使っていない。
-/

namespace ABC3.Found

open Polynomial

/-- **Henselian 局所環では、剰余が `1` の元は `n` 乗根を持つ**
(`n` が単数であるとき)。

`X^n − u` の剰余は `X^n − 1` で、`1` はその根であり、微分の値は `n`
——`n` が `A` の単数なら剰余体で `0` でないので、Hensel が根を持ち上げる。

論文にも我々のモデルにも固有でない**一般の**結果
(`Found/HenselianSplits.lean`・`Found/Teichmuller.lean` と同じ位置づけ)。 -/
theorem exists_pow_eq_of_residue_eq_one {A : Type*} [CommRing A] [HenselianLocalRing A]
    {n : ℕ} (hn : n ≠ 0) (hnu : IsUnit ((n : ℕ) : A)) {u : A}
    (hu : IsLocalRing.residue A u = 1) : ∃ v : A, v ^ n = u := by
  set G : Polynomial A := X ^ n - C u with hG
  have hGm : G.Monic := monic_X_pow_sub_C u hn
  have hGmap : G.map (IsLocalRing.residue A) = X ^ n - 1 := by rw [hG]; simp [hu]
  have hroot : Polynomial.eval (1 : IsLocalRing.ResidueField A)
      (G.map (IsLocalRing.residue A)) = 0 := by rw [hGmap]; simp
  have hder : Polynomial.eval (1 : IsLocalRing.ResidueField A)
      (Polynomial.derivative (G.map (IsLocalRing.residue A))) ≠ 0 := by
    rw [hGmap]
    simp only [derivative_sub, derivative_X_pow, derivative_one, sub_zero, eval_mul, eval_C,
      eval_pow, eval_X, one_pow, mul_one]
    have hc : ((n : ℕ) : IsLocalRing.ResidueField A) = IsLocalRing.residue A ((n : ℕ) : A) := by
      simp
    rw [hc]
    exact (IsLocalRing.residue_ne_zero_iff_isUnit _).mpr hnu
  obtain ⟨v, hv, -⟩ := exists_root_of_residue_root G hGm 1 hroot hder
  refine ⟨v, ?_⟩
  have hz : Polynomial.eval v G = 0 := hv
  rw [hG] at hz
  simp only [eval_sub, eval_pow, eval_X, eval_C, sub_eq_zero] at hz
  exact hz

end ABC3.Found

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Found
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-- **`p` と素な位数の 1 の冪根は `K^ur` の中で `n` 乗根を持つ**(`p ∤ n`)。

`p ∤ n·m` なので `μ_{n·m} ⊆ K^ur`(`exists_isPrimitiveRoot_mem_unramifiedClosure`)。
原始 `n·m` 乗根 `ξ` を取れば `ξ^n` は原始 `m` 乗根なので、`Z^m = 1` なる `Z` は
`(ξ^n)^j` の形に書け、`(ξ^j)^n = Z`。 -/
theorem exists_pow_eq_of_pow_eq_one (K : PAdicLocalField p) {n m : ℕ}
    (hn : ¬ p ∣ n) (hm : ¬ p ∣ m) (hm0 : m ≠ 0) {Z : K.closure} (hZ : Z ^ m = 1) :
    ∃ w : K.closure, w ∈ unramifiedClosure K ∧ w ^ n = Z := by
  have hn0 : n ≠ 0 := by rintro rfl; exact hn (dvd_zero p)
  have hnp : ¬ p ∣ n * m := by
    intro hcon
    rcases (Nat.Prime.dvd_mul (Fact.out : p.Prime)).mp hcon with h | h
    · exact hn h
    · exact hm h
  obtain ⟨ξ, hξmem, hξ⟩ := exists_isPrimitiveRoot_mem_unramifiedClosure K (n * m)
    (Nat.mul_pos (Nat.pos_of_ne_zero hn0) (Nat.pos_of_ne_zero hm0)) hnp
  have hξn : IsPrimitiveRoot (ξ ^ n) m := by
    have h := hξ.pow_of_dvd hn0 (Dvd.intro _ rfl)
    rwa [Nat.mul_div_cancel_left _ (Nat.pos_of_ne_zero hn0)] at h
  haveI : NeZero m := ⟨hm0⟩
  obtain ⟨j, -, hj⟩ := hξn.eq_pow_of_pow_eq_one hZ
  refine ⟨ξ ^ j, pow_mem hξmem j, ?_⟩
  rw [← pow_mul, mul_comm, pow_mul]
  exact hj

/-- 拡大体の整数環の剰余体の位数は `p` で割り切れる——`q_L = q^f = p^{ef}`
(`residueDegree_eq_residueCard_pow` + `residueCard_isPrimePow`、`f ≠ 0` は
`inertiaDegree_ne_zero`)。`p ∤ q_L − 1` を出すために要る。 -/
theorem prime_dvd_card_residueField_adjoinIntegers (K : PAdicLocalField p) (x : K.closure) :
    p ∣ Nat.card (IsLocalRing.ResidueField (adjoinIntegers K x)) := by
  obtain ⟨f, hf, hqf⟩ := residueCard_isPrimePow K
  have hQ : Nat.card (IsLocalRing.ResidueField (adjoinIntegers K x))
      = Nat.card 𝓀[K.carrier] ^ inertiaDegree K x := residueDegree_eq_residueCard_pow K x
  have hq : Nat.card 𝓀[K.carrier] = p ^ f := hqf
  rw [hQ, hq, ← pow_mul]
  exact dvd_pow_self p (by have := inertiaDegree_ne_zero K x; positivity)

/-- **★★★★★★★★★★★★★★★★★★★★`K^ur` のノルム `1` の元は `K^ur` の中に
`n` 乗根を持つ(`p ∤ n`)**——経路 C のノード 1
`units_isNthPow_unramifiedClosure`。

`𝒪_{K^ur}` を構成せず `‖x‖ = 1` で書いてあることについては、この
ファイル冒頭の docstring を参照。

証明は Teichmüller 分解 `x = ζ · u` で「1 の冪根の部分」と
「主単数の部分」に分け、前者は `μ_{n(Q−1)} ⊆ K^ur`(在庫)で、
後者は `adjoinIntegers K x` の中の Hensel で `n` 乗根を作る。 -/
theorem exists_pow_eq_mem_unramifiedClosure (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n)
    {x : K.closure} (hx : x ∈ unramifiedClosure K) (hnorm : ‖x‖ = 1) :
    ∃ y : K.closure, y ∈ unramifiedClosure K ∧ y ^ n = x := by
  have hn0 : n ≠ 0 := by rintro rfl; exact hn (dvd_zero p)
  -- `x ∈ K^ur` は `K(x)/K` 自身が不分岐、ということ
  have hux : IsUnramifiedAdjoin K x := (mem_unramifiedClosure_iff_isUnramified K x).mp hx
  haveI : Fintype (IsLocalRing.ResidueField (adjoinIntegers K x)) := Fintype.ofFinite _
  have hxL : ‖(⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩ :
      IntermediateField.adjoin K.carrier ({x} : Set K.closure))‖ = 1 := hnorm
  set a : adjoinIntegers K x :=
    ⟨⟨x, IntermediateField.mem_adjoin_simple_self K.carrier x⟩, le_of_eq hxL⟩ with hadef
  have ha : IsUnit a := isUnit_of_norm_eq_one K x a hxL
  -- `adjoinIntegers K x → K.closure`(1 層ずつ、`adjoinField` 側へは行かない)
  set ι : adjoinIntegers K x →+* K.closure :=
    (algebraMap (IntermediateField.adjoin K.carrier ({x} : Set K.closure)) K.closure).comp
      (SubringClass.subtype (adjoinIntegers K x)) with hιdef
  have hιa : ι a = x := rfl
  have hιmem : ∀ t : adjoinIntegers K x, ι t ∈ unramifiedClosure K := fun t =>
    adjoin_le_unramifiedClosure K hux
      (t : IntermediateField.adjoin K.carrier ({x} : Set K.closure)).2
  -- 剰余体の位数 `Q`(`p ∣ Q` なので `p ∤ Q − 1`)
  set Q := Fintype.card (IsLocalRing.ResidueField (adjoinIntegers K x)) with hQdef
  have hQ2 : 2 ≤ Q := Fintype.one_lt_card
  have hpQ : p ∣ Q := by
    rw [hQdef, ← Nat.card_eq_fintype_card]
    exact prime_dvd_card_residueField_adjoinIntegers K x
  have hpQ1 : ¬ p ∣ (Q - 1) := by
    intro hcon
    have h1 : p ∣ (Q - (Q - 1)) := Nat.dvd_sub hpQ hcon
    rw [show Q - (Q - 1) = 1 by omega] at h1
    exact (Nat.Prime.one_lt (Fact.out : p.Prime)).ne' (Nat.dvd_one.mp h1)
  -- Teichmüller 分解 `x = ζ · u`(`ζ^{Q−1} = 1`、`u` の剰余は `1`)
  have hb0 : IsLocalRing.residue (adjoinIntegers K x) a ≠ 0 :=
    (IsLocalRing.residue_ne_zero_iff_isUnit a).mpr ha
  obtain ⟨ζ, hζpow, hζres⟩ :=
    exists_teichmullerRep (IsLocalRing.residue (adjoinIntegers K x) a) hb0
  obtain ⟨B, hB⟩ := IsUnit.of_pow_eq_one hζpow (by omega : Q - 1 ≠ 0)
  set u : adjoinIntegers K x := a * (B⁻¹ : (adjoinIntegers K x)ˣ) with hudef
  have hBres : IsLocalRing.residue (adjoinIntegers K x) a
      = IsLocalRing.residue (adjoinIntegers K x) (B : adjoinIntegers K x) := by
    rw [hB, hζres]
  have hures : IsLocalRing.residue (adjoinIntegers K x) u = 1 := by
    rw [hudef, map_mul, hBres, ← map_mul, ← Units.val_mul, mul_inv_cancel, Units.val_one, map_one]
  -- 主単数の部分は `adjoinIntegers K x` の中の Hensel で `n` 乗根を持つ
  obtain ⟨v, hv⟩ := exists_pow_eq_of_residue_eq_one hn0
    (isUnit_natCast_adjoinIntegers K x hn) hures
  have hav : v ^ n * (B : adjoinIntegers K x) = a := by
    rw [hv, hudef, mul_assoc, ← Units.val_mul, inv_mul_cancel, Units.val_one, mul_one]
  -- `K.closure` の側へ移す
  have hmain : (ι v) ^ n * ι (B : adjoinIntegers K x) = x := by
    have h := congrArg ι hav
    rwa [map_mul, map_pow, hιa] at h
  have hZ : (ι (B : adjoinIntegers K x)) ^ (Q - 1) = 1 := by
    rw [← map_pow, show (B : adjoinIntegers K x) ^ (Q - 1) = 1 by rw [hB]; exact hζpow, map_one]
  -- 1 の冪根の部分は `μ_{n(Q−1)} ⊆ K^ur` から `n` 乗根を持つ
  obtain ⟨w, hwmem, hw⟩ := exists_pow_eq_of_pow_eq_one K hn hpQ1 (by omega : Q - 1 ≠ 0) hZ
  refine ⟨w * ι v, mul_mem hwmem (hιmem v), ?_⟩
  rw [mul_pow, hw, mul_comm]
  exact hmain

end ABC3.Found.PGC
