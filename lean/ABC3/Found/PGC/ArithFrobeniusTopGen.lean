import ABC3.Found.PGC.UnramifiedResidueField

/-!
# 算術 Frobenius は整合的な位相的生成元でもある(`sorry` 無し)

経路 Λ の節点 Λ5b′。着地済みの 2 本のあいだに空いていた穴を塞ぐ。

## なぜ要るか(2 本の不整合)

* `UnramifiedZhat.lean::unramifiedClosureGalEquivZHat`(`Gal(K^ur/K) ≃ₜ* Ẑ`)は
  `exists_coherentFrobenius` の `Classical.choose` を使う。そこで取れるのは
  **位相的生成元一般**であり、剰余体上の作用は `z ↦ z^{q^k}`(`k ∈ Ẑ^×`)でよい。
* Λ6(Dwork)は**算術 Frobenius**(`z ↦ z^q` ちょうど)でなければ回らない。
  `UnramifiedResidueField.lean` のモジュール docstring はこの違いを訂正として
  明記し、`exists_arithFrobenius` を別途用意した——ただしそれが
  **位相的生成元でもある**ことは「必要になったら節点を立てること」として
  保留していた(同ファイル docstring の逸脱欄)。
* ⇒ `Ẑ` の同定が `Ẑ^×` のぶん不定なままだと、Λ7 / Λ8 で Artin 写像を組むときに
  正規化が衝突する。本ファイルは**同じ 1 つの元**について両方の性質を示し、
  その穴を塞ぐ。

## 到達点

| | 主張 |
|---|---|
| 主定理 | `exists_arithFrobenius_isCoherent`:剰余体上 `z ↦ z^q` で、かつ**すべての段の生成元**である `σ` が在る |
| 系 | `arithFrobenius` / `arithFrobenius_residue` / `arithFrobenius_mem_unramLevelGeneratorSet` |
| 系 | `zhatMulEquivUnramGalArith`:`Ẑ ≃ₜ* Gal(K^ur/K)` で `1 ↦ 算術 Frobenius`(正規化された同型) |

★`σ` は主定理の `∃` の内側に閉じ込めてある。2 つの性質を**同じ `σ`** について
言うのが要点であり、別々の `∃` に分けたら主定理は無意味になる。

## 筋(位数を剰余体で測る)

`σ` を `exists_arithFrobenius` で取る。第 1 の性質はそのまま。第 2 は
`unramLevelGeneratorSet_eq_preimage` により
「`orderOf (unramLevelHom K N σ) = N`」に帰着する。`Multiplicative (ZMod N)` は
位数 `N` なので `orderOf ∣ N` は自動。逆向き `N ≤ orderOf` が本体で、次のように
剰余体で測る。

1. `σ^k` の剰余体上の作用は `z ↦ z^{q^k}`(`residue_unramGalInt_pow`、`k` の帰納法)。
2. `σ^k` が段 `N`(`unramLevel K N`)を点ごとに固定するとすると、段 `N` の
   整数環の剰余(`𝓀_{K^ur}` の中に `q^N` 元の部分集合を成す——目標 1
   `card_residueField_unramLevel` と `adjResidueToUnram` の単射性)は
   すべて `s^{q^k} = s` を満たす。
3. `X^{q^k} - X` の根はたかだか `q^k` 個(`card_le_of_pow_eq_self`)なので
   `q^N ≤ q^k`、すなわち `N ≤ k`(`le_of_pow_mem_fixingSubgroup_unramLevel`)。
4. `k := orderOf (unramLevelHom K N σ)` に当てる(`k ≠ 0` は群が有限だから)。

★段 `N` の整数の像が `K^ur` の中で本当に `unramLevel K N` に入ることは
`UnramifiedGalCharCount.lean::mem_adjoin_of_val_mem`(`K.closure` の中で
`K⟮z⟯` に入る `K^ur` の元は `K^ur` の中でも `K⟮z⟯` に入る)で渡す。
★ガロア作用と剰余体作用を繋ぐ橋は**在庫にあった**——
`unramGalResidue_residue`(`UnramifiedResidueField.lean`)。本ファイルが新しく
作ったのはその「`k` 乗版」(`residue_unramGalInt_pow`)だけである。

## 在庫でどこまで済んだか(2026-09-06 の実測)

| 要るもの | 在庫 |
|---|---|
| 剰余体上の Frobenius 作用 | `unramGalResidue_residue`(本プロジェクト) |
| 段の剰余体の元の個数 `q^N` | `card_residueField_unramLevel`(本プロジェクト) |
| 段の生成元集合 ↔ 逆像 | `unramLevelGeneratorSet_eq_preimage`(本プロジェクト) |
| `K^ur` の中への降下 | `mem_adjoin_of_val_mem`(本プロジェクト) |
| 根の個数 | `Polynomial.card_roots'` |
| 位数と巡回性 | `orderOf_dvd_natCard` / `Subgroup.eq_top_of_card_eq` / `mem_powers_iff_mem_zpowers` |

★`FiniteField.orderOf_frobeniusAlgEquivOfAlgebraic`(mathlib、有限次の
Frobenius の位数 = `[L:K]`)は**使っていない**。段の剰余体を `𝔽_{q^N}` として
mathlib の有限体と同定する配管(`GaloisField` との同型)を組むより、
`q^N` 元の集合が `X^{q^k} - X` の根に入るという**数え上げ**のほうが安いため。
これは `UnramifiedResidueField.lean` が目標 2 で取ったのと同じ筋である。

## ★設計上の注意(守ったこと)

* **既存ファイルを書き換えていない**(`UnramifiedZhat.lean` /
  `UnramifiedResidueField.lean` は import のみ)。
* **`inertia` を経由していない**。`Interface` の `SubgroupCorrespondence` /
  `ResidueCardinality` は主張にも証明にも現れない。
* **`Abelianization` を使っていない**。
* **結論に自由なパラメータを出していない**——主定理の型は `K` にしか依存せず、
  `σ` は `∃` の内側にある。`arithFrobenius K` はその `Classical.choose`。

## 逸脱(記録)

* 段の剰余体を `𝔽_{q^N}`(`GaloisField p N`)と**同定していない**。
  使うのは「元の個数が `q^N`」だけである(`UnramifiedResidueField.lean` の
  逸脱欄と同じ立場)。下流が「`𝔽_{q^N}` そのもの」を要求したら節点を立てること。
* 算術 Frobenius の**一意性**は主張していない。実際 `Gal(K^ur/K) → Aut(𝓀_{K^ur})`
  は単射なので一意だが、本ファイルはその単射性を要さない
  (下流は `arithFrobenius K` という 1 つの選択を共有すれば足りる)。
  必要になったら節点を立てること。
* `zhatMulEquivUnramGalArith K` は `unramifiedClosureGalEquivZHat K`
  (`UnramifiedZhat.lean`)と**別の同型**である(生成元の取り方が違う)。
  Λ7 / Λ8 で Artin 写像を正規化するときは本ファイルの側を使うこと。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open CategoryTheory ProfiniteGrp ProfiniteGrp.ProfiniteCompletion
open scoped Valued

variable {p : ℕ} [Fact p.Prime]

/-! ## 1. `𝒪_{K^ur}` への制限は群準同型 -/

/-- `Gal(K^ur/K)` の積は `𝒪_{K^ur}` の上で合成。 -/
theorem unramGalInt_mul (K : PAdicLocalField p) (σ τ : unramGal K)
    (w : ↥(unramifiedClosureInt K)) :
    unramGalInt K (σ * τ) w = unramGalInt K σ (unramGalInt K τ w) := by
  apply Subtype.ext
  rw [unramGalInt_coe, unramGalInt_coe, unramGalInt_coe]
  rfl

/-- 単位元は `𝒪_{K^ur}` の上で恒等。 -/
theorem unramGalInt_one (K : PAdicLocalField p) (w : ↥(unramifiedClosureInt K)) :
    unramGalInt K 1 w = w := by
  apply Subtype.ext
  rw [unramGalInt_coe]
  rfl

/-- **算術 Frobenius の `k` 乗は剰余体で `z ↦ z^{q^k}`**。

`k` の帰納法。橋は `unramGalResidue_residue`(在庫)。 -/
theorem residue_unramGalInt_pow (K : PAdicLocalField p) {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) (k : ℕ)
    (w : ↥(unramifiedClosureInt K)) :
    IsLocalRing.residue _ (unramGalInt K (σ ^ k) w)
      = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier] ^ k) := by
  induction k generalizing w with
  | zero => rw [pow_zero, unramGalInt_one, pow_zero, pow_one]
  | succ k ih =>
    rw [pow_succ, unramGalInt_mul, ih, ← unramGalResidue_residue, hσ, ← pow_mul, pow_succ,
      mul_comm]

/-! ## 2. `X^m - X` の根の数え上げ

`UnramifiedResidueField.lean::mem_of_pow_eq_self_of_card_eq` は
「ちょうど `m` 元なら根を尽くす」形。ここで要るのは**上からの評価**だけなので
別に書く(証明は同じ道具)。 -/

open Polynomial in
/-- **`s^m = s` を満たす元の有限集合は `m` 元以下**——`X^m - X` の根だから。

退化の自己検査:`hm`(`1 < m`)を落とすと `m = 1` で `X - X = 0` になり
`card_roots'` が使えない(実際 `m = 1` では任意の `T` が条件を満たすので**偽**)。 -/
theorem card_le_of_pow_eq_self {F : Type*} [Field F] [DecidableEq F] {m : ℕ} (hm : 1 < m)
    {T : Finset F} (hT : ∀ s ∈ T, s ^ m = s) : T.card ≤ m := by
  set P : F[X] := X ^ m - X with hP
  have hdeg : (X : F[X]).natDegree < (X ^ m : F[X]).natDegree := by
    rw [Polynomial.natDegree_X, Polynomial.natDegree_X_pow]; omega
  have hPmonic : P.Monic := by
    refine (Polynomial.monic_X_pow m).sub_of_left ?_
    rw [Polynomial.degree_X, Polynomial.degree_X_pow]
    exact_mod_cast hm
  have hPdeg : P.natDegree = m := by
    rw [hP, Polynomial.natDegree_sub_eq_left_of_natDegree_lt hdeg, Polynomial.natDegree_X_pow]
  have hsub : T ⊆ P.roots.toFinset := by
    intro s hs
    rw [Multiset.mem_toFinset, Polynomial.mem_roots hPmonic.ne_zero, Polynomial.IsRoot.def, hP]
    simp [hT s hs]
  exact le_trans (Finset.card_le_card hsub)
    (le_trans (Multiset.toFinset_card_le _) (le_trans (Polynomial.card_roots' P) (le_of_eq hPdeg)))

/-! ## 3. 段 `N` の整数の像は段 `N` に入る -/

/-- `𝒪_{K(x_N)} → 𝒪_{K^ur}` の像は `K^ur` の中でも `unramLevel K N` に入る。

`mem_adjoin_of_val_mem`(在庫)で `K.closure` 側の所属を `K^ur` 側へ降ろす。 -/
theorem coe_unramIntHomOfLe_mem_unramLevel (K : PAdicLocalField p) {N : ℕ}
    (hle : IntermediateField.adjoin K.carrier ({unramLevelValGen K N} : Set K.closure)
      ≤ unramifiedClosure K)
    (z : adjoinIntegers K (unramLevelValGen K N)) :
    ((unramIntHomOfLe K hle z : ↥(unramifiedClosureInt K)) : ↥(unramifiedClosure K))
      ∈ unramLevel K N := by
  rw [unramLevel]
  refine mem_adjoin_of_val_mem K (unramLevelGen K N) _ ?_
  rw [unramIntHomOfLe_apply]
  exact (z : ↥(IntermediateField.adjoin K.carrier
    ({unramLevelValGen K N} : Set K.closure))).2

/-! ## 4. 位数の下からの評価 -/

/-- **★★★★★★★★★★★★★★算術 Frobenius の `k` 乗が段 `N` を固定するなら `N ≤ k`**。

段 `N` の剰余体は `q^N` 元(`card_residueField_unramLevel`)で、その像は
`𝓀_{K^ur}` の中で `q^N` 元。`σ^k` がそれらを固定すると、どれも `s^{q^k} = s` を
満たすので `q^N ≤ q^k`。

退化の自己検査:`hk0`(`k ≠ 0`)は落とせない——`σ^0 = 1` は任意の段を固定するが
`N ≤ 0` は `hN` に反する。 -/
theorem le_of_pow_mem_fixingSubgroup_unramLevel (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0)
    {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]))
    {k : ℕ} (hk0 : k ≠ 0) (hk : σ ^ k ∈ (unramLevel K N).fixingSubgroup) : N ≤ k := by
  classical
  have hle : IntermediateField.adjoin K.carrier
      ({unramLevelValGen K N} : Set K.closure) ≤ unramifiedClosure K :=
    adjoin_le_unramifiedClosure K (isUnramified_unramLevelValGen K hN)
  haveI : Fintype (IsLocalRing.ResidueField (adjoinIntegers K (unramLevelValGen K N))) :=
    Fintype.ofFinite _
  have hcard : Fintype.card (IsLocalRing.ResidueField (adjoinIntegers K (unramLevelValGen K N)))
      = Nat.card 𝓀[K.carrier] ^ N := by
    rw [← Nat.card_eq_fintype_card]; exact card_residueField_unramLevel K hN
  set φ := adjResidueToUnram K hle with hφ
  set T : Finset (unramResidueField K) := Finset.image φ Finset.univ with hT
  have hTcard : T.card = Nat.card 𝓀[K.carrier] ^ N := by
    rw [hT, Finset.card_image_of_injective _ φ.injective, Finset.card_univ, hcard]
  have hTpow : ∀ s ∈ T, s ^ (Nat.card 𝓀[K.carrier] ^ k) = s := by
    intro s hs
    rw [hT, Finset.mem_image] at hs
    obtain ⟨u, -, rfl⟩ := hs
    obtain ⟨z, rfl⟩ := IsLocalRing.residue_surjective u
    rw [hφ, adjResidueToUnram_residue]
    set w := unramIntHomOfLe K hle z with hw
    have hfix : unramGalInt K (σ ^ k) w = w := by
      apply Subtype.ext
      rw [unramGalInt_coe]
      exact (IntermediateField.mem_fixingSubgroup_iff _ _).mp hk _
        (coe_unramIntHomOfLe_mem_unramLevel K hle z)
    rw [← residue_unramGalInt_pow K hσ k w, hfix]
  have hbound := card_le_of_pow_eq_self (Nat.one_lt_pow hk0 (one_lt_residueCard K)) hTpow
  rw [hTcard] at hbound
  exact (Nat.pow_le_pow_iff_right (one_lt_residueCard K)).mp hbound

/-- **算術 Frobenius の `G/P_N` での位数はちょうど `N`**。

`orderOf ∣ N` は `Multiplicative (ZMod N)` が位数 `N` だから自動、
`N ≤ orderOf` が上の剰余体での評価。 -/
theorem orderOf_unramLevelHom_arithFrobenius (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0)
    {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) :
    orderOf (unramLevelHom K N σ) = N := by
  haveI : NeZero N := ⟨hN⟩
  have hcardG : Nat.card (Multiplicative (ZMod N)) = N := by simp [Nat.card_eq_fintype_card]
  have hdvd : orderOf (unramLevelHom K N σ) ∣ N := by
    have h := orderOf_dvd_natCard (unramLevelHom K N σ)
    rwa [hcardG] at h
  have hpos : orderOf (unramLevelHom K N σ) ≠ 0 := (orderOf_pos _).ne'
  have hmem : σ ^ orderOf (unramLevelHom K N σ) ∈ (unramLevel K N).fixingSubgroup := by
    rw [← (unramLevel_spec K hN).2.2.2, MonoidHom.mem_ker, map_pow, pow_orderOf_eq_one]
  exact le_antisymm (Nat.le_of_dvd (Nat.pos_of_ne_zero hN) hdvd)
    (le_of_pow_mem_fixingSubgroup_unramLevel K hN hσ hpos hmem)

/-- **算術 Frobenius は各段の生成元**。位数が `N` で群が位数 `N` の巡回群だから。 -/
theorem mem_unramLevelGeneratorSet_of_arithFrobenius (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0)
    {σ : unramGal K}
    (hσ : ∀ w : ↥(unramifiedClosureInt K),
      unramGalResidue K σ (IsLocalRing.residue _ w)
        = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) :
    σ ∈ unramLevelGeneratorSet K N := by
  haveI : NeZero N := ⟨hN⟩
  have hcardG : Nat.card (Multiplicative (ZMod N)) = N := by simp [Nat.card_eq_fintype_card]
  have htop : Subgroup.zpowers (unramLevelHom K N σ) = ⊤ := by
    apply Subgroup.eq_top_of_card_eq
    rw [Nat.card_zpowers, orderOf_unramLevelHom_arithFrobenius K hN hσ, hcardG]
  rw [unramLevelGeneratorSet_eq_preimage K hN]
  simp only [Set.mem_preimage, Set.mem_setOf_eq]
  intro b
  have hb : b ∈ Subgroup.zpowers (unramLevelHom K N σ) := htop ▸ Subgroup.mem_top b
  rw [← mem_powers_iff_mem_zpowers, Submonoid.mem_powers_iff] at hb
  exact hb

/-! ## 5. 主定理 -/

/-- **★★★★★★★★★★★★★★★★★★★★★★★★(Λ5b′)算術 Frobenius は整合的な
位相的生成元でもある**——`Gal(K^ur/K)` には、

* `𝓀_{K^ur}` 上で `z ↦ z^q` として作用し(Λ6 = Dwork が要求する側)、
* **同時に**すべての段 `G/P_N` を生成する(Λ5 = `Ẑ` との同定が要求する側)

元 `σ` がある。

**証明**:`exists_arithFrobenius`(`UnramifiedResidueField.lean`)で `σ` を取り、
`mem_unramLevelGeneratorSet_of_arithFrobenius` を当てる。後者は
`orderOf (unramLevelHom K N σ) = N` を剰余体の元の数え上げで出す。

退化の自己検査。

* `z ↦ z^q` を `z ↦ z^{q^k}`(`k ≥ 2`)に替えると第 2 の性質は**偽**——
  段 `k` での位数が `1` になる(`k ∈ Ẑ^×` が要る所以)。
* 2 つの性質を別々の `∃` に分けると、`Λ5` と `Λ6` が**別の元**を指すので
  本定理は無意味になる。★`σ` は `∃` の内側に閉じ込めてある。
* `K^ur` を `K.closure` に替えると第 2 の性質は**偽**(`Γ_K` は非可換)。 -/
theorem exists_arithFrobenius_isCoherent (K : PAdicLocalField p) :
    ∃ σ : unramGal K,
      (∀ w : ↥(unramifiedClosureInt K),
        unramGalResidue K σ (IsLocalRing.residue _ w)
          = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier])) ∧
      (∀ N : ℕ, N ≠ 0 → σ ∈ unramLevelGeneratorSet K N) := by
  obtain ⟨σ, hσ⟩ := exists_arithFrobenius K
  exact ⟨σ, hσ, fun N hN => mem_unramLevelGeneratorSet_of_arithFrobenius K hN hσ⟩

/-! ## 6. 下流のための名前

Λ7 / Λ8 が「同じ元」を指せるように、選択を 1 か所に閉じ込めて名前を付ける。 -/

/-- **算術 Frobenius**(選択)。剰余体上 `z ↦ z^q` で、かつ位相的生成元。 -/
noncomputable def arithFrobenius (K : PAdicLocalField p) : unramGal K :=
  (exists_arithFrobenius_isCoherent K).choose

theorem arithFrobenius_residue (K : PAdicLocalField p) (w : ↥(unramifiedClosureInt K)) :
    unramGalResidue K (arithFrobenius K) (IsLocalRing.residue _ w)
      = (IsLocalRing.residue _ w) ^ (Nat.card 𝓀[K.carrier]) :=
  (exists_arithFrobenius_isCoherent K).choose_spec.1 w

theorem arithFrobenius_mem_unramLevelGeneratorSet (K : PAdicLocalField p) (N : ℕ) (hN : N ≠ 0) :
    arithFrobenius K ∈ unramLevelGeneratorSet K N :=
  (exists_arithFrobenius_isCoherent K).choose_spec.2 N hN

/-- `σ^k` が段 `N` を固定するのは `N ∣ k` のときに限る(算術 Frobenius 版)。 -/
theorem zpow_arithFrobenius_mem_fixingSubgroup_iff (K : PAdicLocalField p) {N : ℕ} (hN : N ≠ 0)
    (k : ℤ) :
    arithFrobenius K ^ k ∈ (unramLevel K N).fixingSubgroup ↔ (N : ℤ) ∣ k :=
  zpow_mem_fixingSubgroup_unramLevel_iff K hN
    (arithFrobenius_mem_unramLevelGeneratorSet K N hN) k

/-- **★★★★★★★★★★★★★★★★★★★★正規化された `Ẑ ≃ₜ* Gal(K^ur/K)`**——
`1 ∈ Ẑ` が算術 Frobenius に対応する。

★`UnramifiedZhat.lean::unramifiedClosureGalEquivZHat` は
`exists_coherentFrobenius` の選択に依存するので `Ẑ^×` のぶん不定だった。
Λ7 / Λ8 で Artin 写像を正規化するときは**こちら**を使うこと。 -/
noncomputable def zhatMulEquivUnramGalArith (K : PAdicLocalField p) : ZHat ≃ₜ* unramGal K :=
  zhatContinuousMulEquivUnramGal K (arithFrobenius_mem_unramLevelGeneratorSet K)

/-- 正規化:`n ∈ ℤ ⊆ Ẑ` は算術 Frobenius の `n` 乗へ行く。 -/
theorem zhatMulEquivUnramGalArith_eta (K : PAdicLocalField p) (n : ℤ) :
    zhatMulEquivUnramGalArith K (ProfiniteCompletion.etaFn _ (Multiplicative.ofAdd n))
      = arithFrobenius K ^ n :=
  ConcreteCategory.congr_hom (ProfiniteCompletion.lift_eta (frobGrpHom K (arithFrobenius K)))
    (Multiplicative.ofAdd n)

end ABC3.Found.PGC
