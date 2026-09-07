import ABC3.Found.PGC.Prop22FixedForm

/-!
# [pGC] Proposition 2.2 —— 残る穴は `𝒪_{K̄}` **だけ**である

`Found/PGC/Prop22FixedForm.lean` は Proposition 2.2(実物に固定した形)を
**ただ 1 つの仮説** `IsometricallyRecoverableClosure p`
(「Prop 2.1 の同変加法同型を等長に取れる」)に還元した。

本ファイルはその穴を**さらに小さくする**:

```
IsometricallyRecoverableClosure p          (等長な K̄ ≃+ K̄′)
  ⟹ UnitBallRecoverableClosure p           (単位球を保つ K̄ ≃+ K̄′)
  ⟺ IntKbarRecoverable (p := p)            (𝒪_{K̄} ≃+ 𝒪_{K̄′} だけ)
  ⟹ IntKbarRecoverable ∧ CompKbarRecoverable   (= Proposition 2.2)
```

すなわち ★★**「`𝒪_{K̄}` が同変加法群として回復できる」ことだけ言えば
Proposition 2.2 は閉じる**。`K̄` 本体も `ℂ_K` も、そこから**無条件に**従う。

## ★★1. なぜこれが進歩か

`Prop22FixedForm.lean` の仮説は「`K̄` の同型を**等長に**取れる」だった。
等長性は `𝒪` の保存より真に強い情報を要求する:

* `𝒪` を保つ加法同型は、`‖x‖ ≤ ‖p‖^k` の球を球へ**正確に**写す
  (`p^k` 倍と可換だから)。しかし `K^al` の値群は `p^ℚ` で**稠密**なので、
  球の対応だけからは `‖φ x‖ = ‖x‖` は出ない —— 出るのは
  `log` の整数部の一致までである。
* 逆に Proposition 2.2 が要求するのは `𝒪_{K̄}` と `ℂ_K` の移送だけで、
  そのどちらにも「正確なノルム」は要らない。

★したがって `IsometricallyRecoverableClosure` は**必要より強い**。
本ファイルはその過剰分を落とす。

★★**原典との関係**: 原文の証明は「有限次拡大 `L/K` に降りて Prop 2.1 を使い、
上付き→下付き→上付きと番号付けを往復して `𝒪_L ⊆ L` を切り出す」。
その**出力**はまさに `𝒪` である。本ファイルの還元は、原文の出力の形と
ぴったり噛み合う —— ★**分岐フィルトレーションが担うべき内容は
`IntKbarRecoverable` の 1 点に集約された**。

## ★★2. 抽象核

* §1 抽象核 C(`transportEquiv`)—— ★**分岐・付値・Galois・ノルムが 1 語も出ない**。
  「加法群 `A` の部分群 `sa` の間の加法同型は、`n` 倍が可逆で `sa` が
  `n` 冪で吸い込むなら、`A` 全体の加法同型に延びる」。
  ★`ψ` が加法的である**だけ**で `ψ (n • w) = n • ψ w` が自動的に従うことが鍵で、
  ここに「`𝒪` の同型は自動的に `p` 倍と可換」という中身が入っている。
* §2 抽象核 D(`continuous_of_unitBall_bound`)—— ★**ノルム付き加法群だけ**。
  「単位球を保つ加法写像は連続」。

## 3. 逸脱の記録

* `Skeleton/PGC/Section2.lean` は書き換えていない。
* `IsometricallyRecoverableClosure` は**消していない** —— `Prop22FixedForm.lean`
  の配線をそのまま残し、本ファイルはその**下**により弱い仮説を挟む。
* `UnitBallRecoverableClosure ⟹ IsometricallyRecoverableClosure` は
  **証明していない**(上で述べた通り、稠密な値群のため一般には出ないと見込む)。
  ★したがって 2 つの仮説の同値性は主張しない。

## ★★4. 気づいた退化 —— `α` に**フィルトレーション両立性が課されていない**

原文の Proposition 2.2 は「`Γ_K` **と** `Γ_K^v`(すべての v > 0)が与えられたとき」と言う。
すなわち回復可能性は「(位相群 + 分岐フィルトレーション)の同型」に対する関手性である。
ところが `Skeleton/PGC/Section2Defs.lean::RecoverableAsAddModule` が量化する `α` は
**位相群の同型だけ**で、`Skeleton/PGC/Section2.lean::prop_2_2` の
`(_RF : RamificationFiltration p)` は**使われていない引数**(先頭の `_`)である。

⇒ ★`IntKbarRecoverable` / `IsometricallyRecoverableClosure` は
**原典より強い主張になっている可能性がある**。Noether の定理により、野性分岐のある
`L/K` では `𝒪_L` は `𝒪_K[Gal(L/K)]`-自由ではないので、Prop 2.1 を通した正規底の議論は
`𝒪` の側にはそのままでは効かない —— 原文が分岐フィルトレーションを**データに入れた**のは
おそらくそのためである。
★**ただしこれが偽であることは示していない。測っていない。**

★★**本ファイルの主結果はこの点に影響されない**: `closureTransport` 以下は
すべて ★**`α` ごと**(`IntKbarTransport α ⟹ …`)の形で述べてある。
後で `α` に「フィルトレーションを保つ」条件を課す形に直しても、そのまま使える。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open UniformSpace

/-! ## §1 抽象核 C —— 部分群の同型を「`n` で割って」全体へ延ばす

★ここには分岐・付値・Galois・ノルム・`PAdicLocalField` は 1 語も出てこない。
出てくるのは「加法群」「その部分群」「自然数 `n` 倍」だけである。 -/

section LocalizationCore

variable {A B : Type*} [AddCommGroup A] [AddCommGroup B] {n : ℕ}
  {SA SB : Type*} [SetLike SA A] [AddSubgroupClass SA A] [SetLike SB B] [AddSubgroupClass SB B]

/-- `n` 倍写像の `k` 回反復は `n^k` 倍。 -/
theorem iterate_nsmulEquiv (mA : A ≃+ A) (hmA : ∀ x, mA x = n • x) (k : ℕ) (x : A) :
    (mA : A → A)^[k] x = (n ^ k) • x := by
  induction k generalizing x with
  | zero => simp
  | succ k ih => rw [Function.iterate_succ_apply, ih, hmA, pow_succ, mul_smul]

/-- 逆向き:`n^k` 倍を `n` 倍の逆で `k` 回割ると戻る。 -/
theorem iterate_symm_nsmulEquiv (mA : A ≃+ A) (hmA : ∀ x, mA x = n • x) (k : ℕ) (x : A) :
    (mA.symm : A → A)^[k] ((n ^ k) • x) = x := by
  rw [← iterate_nsmulEquiv mA hmA k]
  exact (Function.LeftInverse.iterate (g := (mA.symm : A → A)) mA.symm_apply_apply k) x

theorem iterate_addEquiv_add (f : B ≃+ B) (k : ℕ) (b c : B) :
    (f : B → B)^[k] (b + c) = (f : B → B)^[k] b + (f : B → B)^[k] c := by
  induction k generalizing b c with
  | zero => rfl
  | succ k ih =>
      rw [Function.iterate_succ_apply, Function.iterate_succ_apply, Function.iterate_succ_apply,
        map_add, ih]

/-- `n` 倍が同型なら、その逆も作用と可換(作用は加法的だから)。 -/
theorem addEquiv_symm_smul {M' : Type*} [Monoid M'] [DistribMulAction M' B]
    (mB : B ≃+ B) (hmB : ∀ y, mB y = n • y) (g : M') (b : B) :
    mB.symm (g • b) = g • mB.symm b := by
  have h : mB (g • mB.symm b) = g • b := by
    rw [hmB, smul_comm, ← hmB, AddEquiv.apply_symm_apply]
  rw [← h, AddEquiv.symm_apply_apply]

theorem iterate_addEquiv_symm_smul {M' : Type*} [Monoid M'] [DistribMulAction M' B]
    (mB : B ≃+ B) (hmB : ∀ y, mB y = n • y) (g : M') (k : ℕ) (b : B) :
    (mB.symm : B → B)^[k] (g • b) = g • (mB.symm : B → B)^[k] b := by
  induction k generalizing b with
  | zero => rfl
  | succ k ih =>
      rw [Function.iterate_succ_apply, Function.iterate_succ_apply, addEquiv_symm_smul mB hmB, ih]

omit [AddSubgroupClass SA A] in
/-- 部分群への作用の coe が両立していれば、部分群は作用で安定(自動)。 -/
theorem smul_mem_of_coe_smul {M : Type*} [Monoid M] [DistribMulAction M A] {sa : SA}
    [SMul M ↥sa] (hA : ∀ (g : M) (w : ↥sa), ((g • w : ↥sa) : A) = g • (w : A))
    (g : M) {x : A} (hx : x ∈ sa) : g • x ∈ sa := by
  have h := SetLike.coe_mem (g • (⟨x, hx⟩ : ↥sa))
  rwa [hA] at h

theorem pow_nsmul_mem_of_le {k : ℕ} {sa : SA} {x : A} (h : (n ^ k) • x ∈ sa) (m : ℕ) :
    (n ^ (k + m)) • x ∈ sa := by
  rw [add_comm, pow_add, mul_smul]
  exact nsmul_mem h (n ^ m)

theorem pow_nsmul_mem_of_le' {k m : ℕ} (hkm : k ≤ m) {sa : SA} {x : A}
    (h : (n ^ k) • x ∈ sa) : (n ^ m) • x ∈ sa := by
  obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_le hkm
  exact pow_nsmul_mem_of_le h d

/-- 「`x` を `n^k` 倍して `ψ` を当て、`k` 回割り戻す」という候補値。 -/
noncomputable def transportVal (mB : B ≃+ B) {sa : SA} {sb : SB} (ψ : ↥sa ≃+ ↥sb)
    (k : ℕ) (x : A) (h : (n ^ k) • x ∈ sa) : B :=
  (mB.symm : B → B)^[k] ((ψ ⟨(n ^ k) • x, h⟩ : ↥sb) : B)

theorem transportVal_succ (mB : B ≃+ B) (hmB : ∀ y, mB y = n • y) {sa : SA} {sb : SB}
    (ψ : ↥sa ≃+ ↥sb) (k : ℕ) (x : A) (h : (n ^ k) • x ∈ sa) (h' : (n ^ (k + 1)) • x ∈ sa) :
    transportVal mB ψ (k + 1) x h' = transportVal mB ψ k x h := by
  have hx : (⟨(n ^ (k + 1)) • x, h'⟩ : ↥sa) = n • (⟨(n ^ k) • x, h⟩ : ↥sa) := by
    apply Subtype.ext
    push_cast
    rw [pow_succ', mul_smul]
  rw [transportVal, transportVal, hx, map_nsmul]
  have hcoe : ((n • ψ ⟨(n ^ k) • x, h⟩ : ↥sb) : B) = mB ((ψ ⟨(n ^ k) • x, h⟩ : ↥sb) : B) := by
    push_cast
    rw [hmB]
  rw [hcoe, Function.iterate_succ_apply, AddEquiv.symm_apply_apply]

theorem transportVal_add_right (mB : B ≃+ B) (hmB : ∀ y, mB y = n • y) {sa : SA} {sb : SB}
    (ψ : ↥sa ≃+ ↥sb) (k m : ℕ) (x : A) (h : (n ^ k) • x ∈ sa) :
    transportVal mB ψ (k + m) x (pow_nsmul_mem_of_le h m) = transportVal mB ψ k x h := by
  induction m with
  | zero => rfl
  | succ m ih =>
      exact (transportVal_succ mB hmB ψ (k + m) x (pow_nsmul_mem_of_le h m) _).trans ih

/-- ★候補値は `k` の取り方に依らない。 -/
theorem transportVal_eq (mB : B ≃+ B) (hmB : ∀ y, mB y = n • y) {sa : SA} {sb : SB}
    (ψ : ↥sa ≃+ ↥sb) (k l : ℕ) (x : A) (hk : (n ^ k) • x ∈ sa) (hl : (n ^ l) • x ∈ sa) :
    transportVal mB ψ k x hk = transportVal mB ψ l x hl := by
  have key : ∀ (a b : ℕ) (ha : (n ^ a) • x ∈ sa) (hb : (n ^ b) • x ∈ sa), a ≤ b →
      transportVal mB ψ b x hb = transportVal mB ψ a x ha := by
    intro a b ha _ hab
    obtain ⟨m, rfl⟩ := Nat.exists_eq_add_of_le hab
    exact transportVal_add_right mB hmB ψ a m x ha
  rcases le_total k l with h | h
  · exact (key k l hk hl h).symm
  · exact key l k hl hk h

/-- **★★★抽象核 C(写像)** —— 部分群の同型 `ψ` を全体へ延ばした写像。 -/
noncomputable def transportFun (mB : B ≃+ B) {sa : SA} {sb : SB} (ψ : ↥sa ≃+ ↥sb)
    (hcof : ∀ x : A, ∃ k : ℕ, (n ^ k) • x ∈ sa) (x : A) : B :=
  transportVal mB ψ (hcof x).choose x (hcof x).choose_spec

theorem transportFun_eq (mB : B ≃+ B) (hmB : ∀ y, mB y = n • y) {sa : SA} {sb : SB}
    (ψ : ↥sa ≃+ ↥sb) (hcof : ∀ x : A, ∃ k : ℕ, (n ^ k) • x ∈ sa) (k : ℕ) (x : A)
    (h : (n ^ k) • x ∈ sa) : transportFun mB ψ hcof x = transportVal mB ψ k x h :=
  transportVal_eq mB hmB ψ _ k x _ h

/-- ★延長は `sa` の上では `ψ` そのもの。 -/
theorem transportFun_of_mem (mB : B ≃+ B) (hmB : ∀ y, mB y = n • y) {sa : SA} {sb : SB}
    (ψ : ↥sa ≃+ ↥sb) (hcof : ∀ x : A, ∃ k : ℕ, (n ^ k) • x ∈ sa) (x : A) (hx : x ∈ sa) :
    transportFun mB ψ hcof x = ((ψ ⟨x, hx⟩ : ↥sb) : B) := by
  have h0 : (n ^ 0) • x ∈ sa := by simpa using hx
  rw [transportFun_eq mB hmB ψ hcof 0 x h0, transportVal]
  have hx0 : (⟨(n ^ 0) • x, h0⟩ : ↥sa) = ⟨x, hx⟩ := Subtype.ext (by simp)
  rw [hx0]
  rfl

theorem transportFun_add (mB : B ≃+ B) (hmB : ∀ y, mB y = n • y) {sa : SA} {sb : SB}
    (ψ : ↥sa ≃+ ↥sb) (hcof : ∀ x : A, ∃ k : ℕ, (n ^ k) • x ∈ sa) (x y : A) :
    transportFun mB ψ hcof (x + y)
      = transportFun mB ψ hcof x + transportFun mB ψ hcof y := by
  obtain ⟨k, hk⟩ := hcof x
  obtain ⟨l, hl⟩ := hcof y
  have hxm : (n ^ max k l) • x ∈ sa := pow_nsmul_mem_of_le' (le_max_left k l) hk
  have hym : (n ^ max k l) • y ∈ sa := pow_nsmul_mem_of_le' (le_max_right k l) hl
  have hxym : (n ^ max k l) • (x + y) ∈ sa := by rw [smul_add]; exact add_mem hxm hym
  rw [transportFun_eq mB hmB ψ hcof _ x hxm, transportFun_eq mB hmB ψ hcof _ y hym,
    transportFun_eq mB hmB ψ hcof _ (x + y) hxym, transportVal, transportVal, transportVal]
  have hsum : (⟨(n ^ max k l) • (x + y), hxym⟩ : ↥sa)
      = ⟨(n ^ max k l) • x, hxm⟩ + ⟨(n ^ max k l) • y, hym⟩ :=
    Subtype.ext (by push_cast; rw [smul_add])
  rw [hsum, map_add]
  have hcoe : ((ψ ⟨(n ^ max k l) • x, hxm⟩ + ψ ⟨(n ^ max k l) • y, hym⟩ : ↥sb) : B)
      = ((ψ ⟨(n ^ max k l) • x, hxm⟩ : ↥sb) : B)
        + ((ψ ⟨(n ^ max k l) • y, hym⟩ : ↥sb) : B) := by push_cast; rfl
  rw [hcoe, iterate_addEquiv_add]

/-- ★`n^k` 倍すると `ψ` の像に戻る。 -/
theorem nsmul_transportFun (mB : B ≃+ B) (hmB : ∀ y, mB y = n • y) {sa : SA} {sb : SB}
    (ψ : ↥sa ≃+ ↥sb) (hcof : ∀ x : A, ∃ k : ℕ, (n ^ k) • x ∈ sa) (k : ℕ) (x : A)
    (h : (n ^ k) • x ∈ sa) :
    (n ^ k) • transportFun mB ψ hcof x = ((ψ ⟨(n ^ k) • x, h⟩ : ↥sb) : B) := by
  rw [transportFun_eq mB hmB ψ hcof k x h, transportVal,
    ← iterate_nsmulEquiv (A := B) mB hmB k]
  exact (Function.LeftInverse.iterate (g := (mB : B → B)) mB.apply_symm_apply k) _

theorem transportFun_leftInverse (mA : A ≃+ A) (mB : B ≃+ B) (hmA : ∀ x, mA x = n • x)
    (hmB : ∀ y, mB y = n • y) {sa : SA} {sb : SB} (ψ : ↥sa ≃+ ↥sb)
    (hcofA : ∀ x : A, ∃ k : ℕ, (n ^ k) • x ∈ sa)
    (hcofB : ∀ y : B, ∃ k : ℕ, (n ^ k) • y ∈ sb) (x : A) :
    transportFun mA ψ.symm hcofB (transportFun mB ψ hcofA x) = x := by
  obtain ⟨k, hk⟩ := hcofA x
  have hmem : (n ^ k) • transportFun mB ψ hcofA x ∈ sb := by
    rw [nsmul_transportFun mB hmB ψ hcofA k x hk]
    exact SetLike.coe_mem _
  rw [transportFun_eq mA hmA ψ.symm hcofB k _ hmem, transportVal]
  have hval : (⟨(n ^ k) • transportFun mB ψ hcofA x, hmem⟩ : ↥sb) = ψ ⟨(n ^ k) • x, hk⟩ :=
    Subtype.ext (nsmul_transportFun mB hmB ψ hcofA k x hk)
  rw [hval, AddEquiv.symm_apply_apply]
  exact iterate_symm_nsmulEquiv mA hmA k x

/-- **★★★抽象核 C** —— 部分群の間の加法同型は全体の加法同型に延びる。 -/
noncomputable def transportEquiv (mA : A ≃+ A) (mB : B ≃+ B) (hmA : ∀ x, mA x = n • x)
    (hmB : ∀ y, mB y = n • y) {sa : SA} {sb : SB} (ψ : ↥sa ≃+ ↥sb)
    (hcofA : ∀ x : A, ∃ k : ℕ, (n ^ k) • x ∈ sa)
    (hcofB : ∀ y : B, ∃ k : ℕ, (n ^ k) • y ∈ sb) : A ≃+ B where
  toFun := transportFun mB ψ hcofA
  invFun := transportFun mA ψ.symm hcofB
  left_inv := transportFun_leftInverse mA mB hmA hmB ψ hcofA hcofB
  right_inv y := by
    have h := transportFun_leftInverse mB mA hmB hmA ψ.symm hcofB hcofA y
    rwa [AddEquiv.symm_symm] at h
  map_add' := transportFun_add mB hmB ψ hcofA

@[simp] theorem transportEquiv_apply (mA : A ≃+ A) (mB : B ≃+ B) (hmA : ∀ x, mA x = n • x)
    (hmB : ∀ y, mB y = n • y) {sa : SA} {sb : SB} (ψ : ↥sa ≃+ ↥sb)
    (hcofA : ∀ x : A, ∃ k : ℕ, (n ^ k) • x ∈ sa)
    (hcofB : ∀ y : B, ∃ k : ℕ, (n ^ k) • y ∈ sb) (x : A) :
    transportEquiv mA mB hmA hmB ψ hcofA hcofB x = transportFun mB ψ hcofA x := rfl

@[simp] theorem transportEquiv_symm_apply (mA : A ≃+ A) (mB : B ≃+ B) (hmA : ∀ x, mA x = n • x)
    (hmB : ∀ y, mB y = n • y) {sa : SA} {sb : SB} (ψ : ↥sa ≃+ ↥sb)
    (hcofA : ∀ x : A, ∃ k : ℕ, (n ^ k) • x ∈ sa)
    (hcofB : ∀ y : B, ∃ k : ℕ, (n ^ k) • y ∈ sb) (y : B) :
    (transportEquiv mA mB hmA hmB ψ hcofA hcofB).symm y = transportFun mA ψ.symm hcofB y := rfl

/-- ★★延長は部分群を**ちょうど**部分群に写す。 -/
theorem mem_iff_transportEquiv (mA : A ≃+ A) (mB : B ≃+ B) (hmA : ∀ x, mA x = n • x)
    (hmB : ∀ y, mB y = n • y) {sa : SA} {sb : SB} (ψ : ↥sa ≃+ ↥sb)
    (hcofA : ∀ x : A, ∃ k : ℕ, (n ^ k) • x ∈ sa)
    (hcofB : ∀ y : B, ∃ k : ℕ, (n ^ k) • y ∈ sb) (x : A) :
    x ∈ sa ↔ transportEquiv mA mB hmA hmB ψ hcofA hcofB x ∈ sb := by
  constructor
  · intro hx
    rw [transportEquiv_apply, transportFun_of_mem mB hmB ψ hcofA x hx]
    exact SetLike.coe_mem _
  · intro hx
    rw [transportEquiv_apply] at hx
    have hback := transportFun_leftInverse mA mB hmA hmB ψ hcofA hcofB x
    rw [transportFun_of_mem mA hmA ψ.symm hcofB _ hx] at hback
    rw [← hback]
    exact SetLike.coe_mem _

/-- ★★延長も同変。 -/
theorem transportFun_smul {M M' : Type*} [Monoid M] [Monoid M'] [DistribMulAction M A]
    [DistribMulAction M' B] (mB : B ≃+ B) (hmB : ∀ y, mB y = n • y) {sa : SA} {sb : SB}
    [SMul M ↥sa] [SMul M' ↥sb] (ψ : ↥sa ≃+ ↥sb)
    (hcof : ∀ x : A, ∃ k : ℕ, (n ^ k) • x ∈ sa)
    (hA : ∀ (g : M) (w : ↥sa), ((g • w : ↥sa) : A) = g • (w : A))
    (hB : ∀ (g' : M') (w : ↥sb), ((g' • w : ↥sb) : B) = g' • (w : B))
    (α : M → M') (heq : ∀ (g : M) (w : ↥sa), ψ (g • w) = α g • ψ w) (g : M) (x : A) :
    transportFun mB ψ hcof (g • x) = α g • transportFun mB ψ hcof x := by
  obtain ⟨k, hk⟩ := hcof x
  have hgk : (n ^ k) • (g • x) ∈ sa := by
    rw [smul_comm]
    exact smul_mem_of_coe_smul hA g hk
  rw [transportFun_eq mB hmB ψ hcof k (g • x) hgk, transportFun_eq mB hmB ψ hcof k x hk,
    transportVal, transportVal]
  have hsub : (⟨(n ^ k) • (g • x), hgk⟩ : ↥sa) = g • ⟨(n ^ k) • x, hk⟩ :=
    Subtype.ext (by rw [hA]; exact smul_comm _ _ _)
  rw [hsub, heq, hB, iterate_addEquiv_symm_smul mB hmB]

end LocalizationCore

/-! ## §2 抽象核 D —— 単位球を保つ加法写像は連続

★ここにも分岐・付値・Galois は 1 語も出てこない。「ノルム付き加法群」だけである。 -/

section UnitBallCore

variable {A B : Type*} [SeminormedAddCommGroup A] [SeminormedAddCommGroup B]

/-- 加法写像の連続性は「0 の近くでの評価」だけで出る。 -/
theorem continuous_of_ball (φ : A →+ B)
    (h : ∀ ε : ℝ, 0 < ε → ∃ δ > 0, ∀ a : A, ‖a‖ < δ → ‖φ a‖ < ε) : Continuous φ := by
  rw [Metric.continuous_iff]
  intro b ε hε
  obtain ⟨δ, hδ, hmain⟩ := h ε hε
  refine ⟨δ, hδ, fun a ha => ?_⟩
  have hd : dist (φ a) (φ b) = ‖φ (a - b)‖ := by rw [map_sub]; exact dist_eq_norm _ _
  rw [hd]
  exact hmain _ (by rwa [← dist_eq_norm])

/-- **★★★抽象核 D(評価)** —— 単位球を保つ加法写像は、半径 `cA^k` の球を
半径 `cB^k` の球へ写す。

`n` 倍が `A` で `cA` 倍、`B` で `cB` 倍にノルムを変え、`A` で `n` 倍が全射なら。
★加法写像は `φ (n • a) = n • φ a` を**自動的に**満たす —— そこが効いている。 -/
theorem norm_le_pow_of_unitBall {n : ℕ} (φ : A →+ B) (hball : ∀ a : A, ‖a‖ ≤ 1 → ‖φ a‖ ≤ 1)
    {cA cB : ℝ} (hcA : 0 < cA) (hcB : 0 ≤ cB) (hnA : ∀ a : A, ‖n • a‖ = cA * ‖a‖)
    (hnB : ∀ b : B, ‖n • b‖ = cB * ‖b‖) (hdiv : ∀ a : A, ∃ a', n • a' = a)
    (k : ℕ) (a : A) (h : ‖a‖ ≤ cA ^ k) : ‖φ a‖ ≤ cB ^ k := by
  induction k generalizing a with
  | zero => simpa using hball a (by simpa using h)
  | succ k ih =>
      obtain ⟨a', rfl⟩ := hdiv a
      have h' : ‖a'‖ ≤ cA ^ k := by
        rw [hnA, pow_succ'] at h
        exact le_of_mul_le_mul_left h hcA
      have hstep := ih a' h'
      rw [map_nsmul, hnB, pow_succ']
      exact mul_le_mul_of_nonneg_left hstep hcB

/-- **★★★抽象核 D(連続性)**。 -/
theorem continuous_of_unitBall_bound {n : ℕ} (φ : A →+ B)
    (hball : ∀ a : A, ‖a‖ ≤ 1 → ‖φ a‖ ≤ 1) {cA cB : ℝ} (hcA : 0 < cA) (hcB : 0 ≤ cB)
    (hcB1 : cB < 1) (hnA : ∀ a : A, ‖n • a‖ = cA * ‖a‖) (hnB : ∀ b : B, ‖n • b‖ = cB * ‖b‖)
    (hdiv : ∀ a : A, ∃ a', n • a' = a) : Continuous φ := by
  refine continuous_of_ball φ (fun ε hε => ?_)
  obtain ⟨k, hk⟩ := exists_pow_lt_of_lt_one hε hcB1
  exact ⟨cA ^ k, pow_pos hcA k, fun a ha =>
    lt_of_le_of_lt (norm_le_pow_of_unitBall φ hball hcA hcB hnA hnB hdiv k a ha.le) hk⟩

end UnitBallCore

/-! ## §3 具体層 —— `K^al` と `𝒪_{K̄}` への代入 -/

variable {p : ℕ} [Fact p.Prime]

/-- 体の非零元による左乗法(加法同型)。★mathlib には `Equiv.mulLeft₀` はあるが
`AddEquiv` 版は無い(実測: `grep -n "mulLeft₀" .cache/mathlib-index.txt` は
`Equiv` / `Homeomorph` / `MeasurableEquiv` / `OrderIso` の 4 つだけを返す)。 -/
noncomputable def mulLeftAddEquivOfNeZero {F : Type*} [Field F] (a : F) (ha : a ≠ 0) : F ≃+ F where
  toFun x := a * x
  invFun x := a⁻¹ * x
  left_inv x := by field_simp
  right_inv x := by field_simp
  map_add' x y := mul_add _ _ _

/-- `K^al` 上の「`p` 倍」加法同型。 -/
noncomputable def closureMulP (K : PAdicLocalField p) : K.closure ≃+ K.closure :=
  mulLeftAddEquivOfNeZero ((p : ℕ) : K.closure) (natCast_p_closure_ne_zero K)

theorem nsmul_p_closure_eq (K : PAdicLocalField p) (x : K.closure) :
    ((p : ℕ) • x : K.closure) = ((p : ℕ) : K.closure) * x := nsmul_eq_mul _ _

theorem closureMulP_apply (K : PAdicLocalField p) (x : K.closure) :
    closureMulP K x = (p : ℕ) • x := (nsmul_p_closure_eq K x).symm

theorem norm_nsmul_p_closure (K : PAdicLocalField p) (x : K.closure) :
    ‖((p : ℕ) • x : K.closure)‖ = ((p : ℕ) : ℝ)⁻¹ * ‖x‖ := by
  rw [nsmul_p_closure_eq, norm_mul, norm_natCast_p_closure]

theorem exists_nsmul_p_closure (K : PAdicLocalField p) (x : K.closure) :
    ∃ y : K.closure, ((p : ℕ) • y : K.closure) = x := by
  have hp : ((p : ℕ) : K.closure) ≠ 0 := natCast_p_closure_ne_zero K
  refine ⟨x / ((p : ℕ) : K.closure), ?_⟩
  rw [nsmul_p_closure_eq]
  field_simp

/-- ★`𝒪_{K̄}` は `p` 冪で `K^al` を吸い込む(抽象核 C の共終性仮説)。 -/
theorem exists_pow_nsmul_mem_absClosureInt (K : PAdicLocalField p) (x : K.closure) :
    ∃ k : ℕ, ((((p : ℕ)) ^ k) • x : K.closure) ∈ absClosureInt K := by
  by_cases hx : x = 0
  · exact ⟨0, by simp [hx]⟩
  have hx0 : 0 < ‖x‖ := norm_pos_iff.mpr hx
  obtain ⟨k, hk⟩ := exists_pow_lt_of_lt_one (inv_pos.mpr hx0) (norm_natCast_p_closure_lt_one K)
  refine ⟨k, ?_⟩
  rw [mem_absClosureInt]
  have hsm : ((((p : ℕ)) ^ k) • x : K.closure) = ((p : ℕ) : K.closure) ^ k * x := by
    rw [nsmul_eq_mul, Nat.cast_pow]
  rw [hsm, norm_mul, norm_pow]
  calc ‖((p : ℕ) : K.closure)‖ ^ k * ‖x‖
      ≤ ‖x‖⁻¹ * ‖x‖ := mul_le_mul_of_nonneg_right hk.le (norm_nonneg _)
    _ = 1 := inv_mul_cancel₀ (ne_of_gt hx0)

/-! ### `𝒪_{K̄}` の移送を `K^al` の移送へ延ばす -/

/-- **α ごとの `𝒪_{K̄}` の移送**という仮説(`IntKbarRecoverable` の α ごとの形)。 -/
def IntKbarTransport {K K' : PAdicLocalField p}
    (α : ContinuousMulEquiv K.absGal K'.absGal) : Prop :=
  ∃ φ : IntKbar K ≃+ IntKbar K', ∀ (g : K.absGal) (x : IntKbar K),
    φ (g • x) = (α.toMulEquiv g) • φ x

/-- ★★★**`𝒪_{K̄} ≃+ 𝒪_{K̄′}` は `K^al ≃+ K̄′` に延びる**(抽象核 C の代入)。 -/
noncomputable def closureTransport {K K' : PAdicLocalField p}
    (ψ : IntKbar K ≃+ IntKbar K') : K.closure ≃+ K'.closure :=
  transportEquiv (closureMulP K) (closureMulP K') (closureMulP_apply K) (closureMulP_apply K')
    ψ (exists_pow_nsmul_mem_absClosureInt K) (exists_pow_nsmul_mem_absClosureInt K')

/-- ★延長は `𝒪_{K̄}` をちょうど `𝒪_{K̄′}` に写す。 -/
theorem mem_iff_closureTransport {K K' : PAdicLocalField p} (ψ : IntKbar K ≃+ IntKbar K')
    (x : K.closure) : x ∈ absClosureInt K ↔ closureTransport ψ x ∈ absClosureInt K' :=
  mem_iff_transportEquiv _ _ _ _ ψ _ _ x

theorem norm_le_one_closureTransport {K K' : PAdicLocalField p} (ψ : IntKbar K ≃+ IntKbar K')
    (x : K.closure) (hx : ‖x‖ ≤ 1) : ‖closureTransport ψ x‖ ≤ 1 := by
  rw [← mem_absClosureInt K']
  rw [← mem_absClosureInt K] at hx
  exact (mem_iff_closureTransport ψ x).mp hx

/-- ★延長は連続(抽象核 D の代入)。★等長性は**要らない**。 -/
theorem continuous_closureTransport {K K' : PAdicLocalField p} (ψ : IntKbar K ≃+ IntKbar K') :
    Continuous (closureTransport ψ) := by
  have hp1 : (1 : ℝ) < ((p : ℕ) : ℝ) := by exact_mod_cast (Fact.out : p.Prime).one_lt
  have hpos : (0 : ℝ) < ((p : ℕ) : ℝ)⁻¹ := inv_pos.mpr (lt_trans zero_lt_one hp1)
  exact continuous_of_unitBall_bound (n := (p : ℕ)) (closureTransport ψ).toAddMonoidHom
    (norm_le_one_closureTransport ψ) hpos hpos.le (inv_lt_one_of_one_lt₀ hp1)
    (norm_nsmul_p_closure K) (norm_nsmul_p_closure K') (exists_nsmul_p_closure K)

theorem continuous_closureTransport_symm {K K' : PAdicLocalField p}
    (ψ : IntKbar K ≃+ IntKbar K') : Continuous (closureTransport ψ).symm :=
  continuous_closureTransport ψ.symm

/-- ★延長も同変(抽象核 C の代入)。 -/
theorem closureTransport_smul {K K' : PAdicLocalField p}
    {α : ContinuousMulEquiv K.absGal K'.absGal} (ψ : IntKbar K ≃+ IntKbar K')
    (hψ : ∀ (g : K.absGal) (w : IntKbar K), ψ (g • w) = (α.toMulEquiv g) • ψ w)
    (g : K.absGal) (x : K.closure) :
    closureTransport ψ (g • x) = (α.toMulEquiv g) • closureTransport ψ x :=
  transportFun_smul (closureMulP K') (closureMulP_apply K') ψ
    (exists_pow_nsmul_mem_absClosureInt K)
    (fun g w => by rw [coe_smul_absClosureInt, smul_closure_def])
    (fun g w => by rw [coe_smul_absClosureInt, smul_closure_def])
    (fun g => α.toMulEquiv g) hψ g x

/-! ### 結論 —— `𝒪_{K̄}` だけで Proposition 2.2 が閉じる -/

/-- ★★**`𝒪_{K̄}` の移送から `K^al` の移送**(= Proposition 2.1)。 -/
theorem closure_transport_of_intKbarTransport {K K' : PAdicLocalField p}
    {α : ContinuousMulEquiv K.absGal K'.absGal} (h : IntKbarTransport α) :
    ∃ φ : K.closure ≃+ K'.closure, ∀ (g : K.absGal) (x : K.closure),
      φ (g • x) = (α.toMulEquiv g) • φ x := by
  obtain ⟨ψ, hψ⟩ := h
  exact ⟨closureTransport ψ, closureTransport_smul ψ hψ⟩

/-- ★★★**`𝒪_{K̄}` の移送から `ℂ_K` の移送**。 -/
theorem compKbar_transport_of_intKbarTransport {K K' : PAdicLocalField p}
    {α : ContinuousMulEquiv K.absGal K'.absGal} (h : IntKbarTransport α) :
    ∃ φ : CompKbar K ≃+ CompKbar K', ∀ (g : K.absGal) (x : CompKbar K),
      φ (g • x) = (α.toMulEquiv g) • φ x := by
  obtain ⟨ψ, hψ⟩ := h
  exact ⟨addEquivCompletion (closureTransport ψ) (continuous_closureTransport ψ)
      (continuous_closureTransport_symm ψ),
    addEquivCompletion_smul _ _ _ (fun g => α.toMulEquiv g) (closureTransport_smul ψ hψ)⟩

/-- ★★★**[pGC] Proposition 2.2(実物に固定した形)—— 仮説は `𝒪_{K̄}` だけでよい**。

`Prop22FixedForm.lean::prop_2_2_real_of_isometric` は「等長な `K̄` の同型」を要求したが、
実際には ★**`𝒪_{K̄}` の同変加法同型があれば十分**である。 -/
theorem prop_2_2_real_of_intKbar (h : IntKbarRecoverable (p := p)) :
    IntKbarRecoverable (p := p) ∧ CompKbarRecoverable (p := p) :=
  ⟨h, fun α => compKbar_transport_of_intKbarTransport (h α)⟩

def prop_2_2_real_of_intKbar.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★**Proposition 2.1 は `𝒪_{K̄}` の回復から従う**。 -/
theorem recoverableAsAddModule_closure_of_intKbar (h : IntKbarRecoverable (p := p)) :
    RecoverableAsAddModule (p := p) (fun K => K.closure) :=
  fun α => closure_transport_of_intKbarTransport (h α)

/-! ## §4 単位球版 —— 等長版と `𝒪` 版のあいだ -/

/-- **α ごとの「単位球版 Proposition 2.1」**。等長性を「単位球の対応」に弱めたもの。 -/
def UnitBallTransport {K K' : PAdicLocalField p}
    (α : ContinuousMulEquiv K.absGal K'.absGal) : Prop :=
  ∃ φ : K.closure ≃+ K'.closure, (∀ x, ‖x‖ ≤ 1 ↔ ‖φ x‖ ≤ 1) ∧
    ∀ (g : K.absGal) (x : K.closure), φ (g • x) = (α.toMulEquiv g) • φ x

/-- **単位球版 Proposition 2.1**(すべての `α` について)。 -/
def UnitBallRecoverableClosure (p : ℕ) [Fact p.Prime] : Prop :=
  ∀ {K K' : PAdicLocalField p} (α : ContinuousMulEquiv K.absGal K'.absGal), UnitBallTransport α

/-- ★等長版 ⟹ 単位球版(真に弱める向き)。 -/
theorem unitBallTransport_of_isometricTransport {K K' : PAdicLocalField p}
    {α : ContinuousMulEquiv K.absGal K'.absGal} (h : IsometricTransport α) :
    UnitBallTransport α := by
  obtain ⟨φ, hn, he⟩ := h
  exact ⟨φ, fun x => by rw [hn], he⟩

/-- ★単位球版 ⟹ `𝒪` 版。 -/
theorem intKbarTransport_of_unitBallTransport {K K' : PAdicLocalField p}
    {α : ContinuousMulEquiv K.absGal K'.absGal} (h : UnitBallTransport α) :
    IntKbarTransport α := by
  obtain ⟨ψ, hball, hequiv⟩ := h
  have hmem : ∀ x : K.closure, x ∈ absClosureInt K ↔ ψ x ∈ absClosureInt K' := by
    intro x; rw [mem_absClosureInt, mem_absClosureInt, hball]
  refine ⟨addEquivOfMemIff ψ (absClosureInt K) (absClosureInt K') hmem, fun g x => ?_⟩
  exact Subtype.ext (by
    rw [coe_addEquivOfMemIff, coe_smul_absClosureInt, coe_smul_absClosureInt,
      coe_addEquivOfMemIff, ← smul_closure_def, hequiv, smul_closure_def])

/-- ★`𝒪` 版 ⟹ 単位球版(逆向き。抽象核 C の延長がそのまま単位球版になる)。 -/
theorem unitBallTransport_of_intKbarTransport {K K' : PAdicLocalField p}
    {α : ContinuousMulEquiv K.absGal K'.absGal} (h : IntKbarTransport α) :
    UnitBallTransport α := by
  obtain ⟨ψ, hψ⟩ := h
  refine ⟨closureTransport ψ, fun x => ?_, closureTransport_smul ψ hψ⟩
  rw [← mem_absClosureInt K, ← mem_absClosureInt K']
  exact mem_iff_closureTransport ψ x

/-- ★★**単位球版と `𝒪` 版は同値**。 -/
theorem unitBallTransport_iff_intKbarTransport {K K' : PAdicLocalField p}
    {α : ContinuousMulEquiv K.absGal K'.absGal} :
    UnitBallTransport α ↔ IntKbarTransport α :=
  ⟨intKbarTransport_of_unitBallTransport, unitBallTransport_of_intKbarTransport⟩

theorem unitBallRecoverableClosure_of_isometric (h : IsometricallyRecoverableClosure p) :
    UnitBallRecoverableClosure p := fun α => unitBallTransport_of_isometricTransport (h α)

/-- ★★**[pGC] Proposition 2.2 —— 単位球版からの還元**。 -/
theorem prop_2_2_real_of_unitBall (h : UnitBallRecoverableClosure p) :
    IntKbarRecoverable (p := p) ∧ CompKbarRecoverable (p := p) :=
  prop_2_2_real_of_intKbar (fun α => intKbarTransport_of_unitBallTransport (h α))

def prop_2_2_real_of_unitBall.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 5, item := "Proposition 2.2", sectionId := "prop-2-2" }

/-- ★★**弱めた鎖が旧結論を再現することの機械検査**。

`Prop22FixedForm.lean::prop_2_2_real_of_isometric` と同じ結論が、
`IsometricallyRecoverableClosure ⟹ UnitBall ⟹ 𝒪` の経路でも出る。
★すなわち本ファイルの仮説は旧仮説より**弱い**(強くはなっていない)。 -/
theorem prop_2_2_real_of_isometric_via_intKbar (h : IsometricallyRecoverableClosure p) :
    IntKbarRecoverable (p := p) ∧ CompKbarRecoverable (p := p) :=
  prop_2_2_real_of_unitBall (unitBallRecoverableClosure_of_isometric h)

/-! ## §5 仮説たちが群体(groupoid)であること

★`IsometricTransport` / `IntKbarTransport` が成り立つ `α` の全体は、
恒等・逆・合成で閉じている。`Prop22FixedForm.lean` の非空虚性
(内部自己同型・体の同型から来る `α`)は、これで**生成された部分群体**の上で成り立つ。 -/

theorem isometricTransport_symm {K K' : PAdicLocalField p}
    {α : ContinuousMulEquiv K.absGal K'.absGal} (h : IsometricTransport α) :
    IsometricTransport α.symm := by
  obtain ⟨φ, hn, he⟩ := h
  refine ⟨φ.symm, fun y => by rw [← hn (φ.symm y), φ.apply_symm_apply], fun g' y => ?_⟩
  have hkey := he (α.symm.toMulEquiv g') (φ.symm y)
  rw [φ.apply_symm_apply] at hkey
  have hgg : α.toMulEquiv (α.symm.toMulEquiv g') = g' := α.apply_symm_apply g'
  rw [hgg] at hkey
  rw [← hkey, φ.symm_apply_apply]

theorem intKbarTransport_refl (K : PAdicLocalField p) :
    IntKbarTransport (ContinuousMulEquiv.refl K.absGal) :=
  ⟨AddEquiv.refl _, fun _ _ => rfl⟩

theorem intKbarTransport_symm {K K' : PAdicLocalField p}
    {α : ContinuousMulEquiv K.absGal K'.absGal} (h : IntKbarTransport α) :
    IntKbarTransport α.symm := by
  obtain ⟨φ, he⟩ := h
  refine ⟨φ.symm, fun g' y => ?_⟩
  have hkey := he (α.symm.toMulEquiv g') (φ.symm y)
  have hgg : α.toMulEquiv (α.symm.toMulEquiv g') = g' := α.apply_symm_apply g'
  rw [φ.apply_symm_apply, hgg] at hkey
  rw [← hkey, φ.symm_apply_apply]

theorem intKbarTransport_trans {K K' K'' : PAdicLocalField p}
    {α : ContinuousMulEquiv K.absGal K'.absGal}
    {β : ContinuousMulEquiv K'.absGal K''.absGal}
    (hα : IntKbarTransport α) (hβ : IntKbarTransport β) :
    IntKbarTransport (α.trans β) := by
  obtain ⟨φ, heφ⟩ := hα
  obtain ⟨χ, heχ⟩ := hβ
  refine ⟨φ.trans χ, fun g x => ?_⟩
  rw [AddEquiv.trans_apply, AddEquiv.trans_apply, heφ, heχ]
  rfl

/-! ### 新しい仮説の非空虚性

★`Prop22FixedForm.lean` §4 の 2 つの場合は、そのまま `IntKbarTransport` の
非空虚性を与える —— ★**新しい(より弱い)仮説も空虚ではない**。 -/

theorem intKbarTransport_of_isometricTransport {K K' : PAdicLocalField p}
    {α : ContinuousMulEquiv K.absGal K'.absGal} (h : IsometricTransport α) :
    IntKbarTransport α := intKbar_transport_of_isometricTransport h

theorem intKbarTransport_inner (K : PAdicLocalField p) (c : K.absGal) :
    IntKbarTransport (innerAbsGalEquiv K c) :=
  intKbarTransport_of_isometricTransport (isometricTransport_inner K c)

/-- ★体の同型から来る `α`(= 自由版の反例が使う α そのもの)では**無条件に**成り立つ。 -/
theorem intKbarTransport_galContinuousMulEquiv {K K' : PAdicLocalField p}
    (β : K.carrier ≃ₐ[ℚ_[p]] K'.carrier) :
    IntKbarTransport (galContinuousMulEquiv β) :=
  intKbarTransport_of_isometricTransport (isometricTransport_galContinuousMulEquiv β)

theorem isometricTransport_trans {K K' K'' : PAdicLocalField p}
    {α : ContinuousMulEquiv K.absGal K'.absGal}
    {β : ContinuousMulEquiv K'.absGal K''.absGal}
    (hα : IsometricTransport α) (hβ : IsometricTransport β) :
    IsometricTransport (α.trans β) := by
  obtain ⟨φ, hnφ, heφ⟩ := hα
  obtain ⟨χ, hnχ, heχ⟩ := hβ
  refine ⟨φ.trans χ, fun x => by rw [AddEquiv.trans_apply, hnχ, hnφ], fun g x => ?_⟩
  rw [AddEquiv.trans_apply, AddEquiv.trans_apply, heφ, heχ]
  rfl

end ABC3.Found.PGC
