/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.TensorPowTriv
import ABC3.Found.Arakelov.TensorMetric
import ABC3.Meta.Claim

/-!
# ★★★★★★★★重なりの一致判定 —— 遷移単元の `N` 乗（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★★★これは何か —— 段 E3a-3 の判定基準

段 E3a-3（大域化）では、チャート `V` ごとに得た**関数** `a ∈ Γ(X, V)` を
`M^{⊗N}` の切断 `secOfFun` に戻して貼り合わせる。
★★**「重なりで一致するか」を関数の言葉で判定する**のが本ファイルである:

    `secOfFun e N a` と `secOfFun e' N a'` が `V ⊓ V'` で一致
      ⟺  `a|_{V⊓V'} · u^N = a'|_{V⊓V'}`   （`u` は `e → e'` の遷移単元）

★★★すなわち**遷移単元がちょうど `N` 乗で効く**。これが「`n` 冪へ移る」議論の
幾何的な意味である——`M` の遷移函数を `N` 乗したものが `M^{⊗N}` の遷移函数だからである。

## ★★★★機構 —— 在庫 3 本＋新規 3 本

| 段 | 道具 |
|---|---|
| 遷移単元の `N` 乗 | `transUnit_tensorTriv`（在庫、`TensorMetric.lean`）を再帰 |
| 制限との可換 | `trivialOfLe` は `tensorTriv` と可換（`rfl`）→ 再帰 |
| 座標の移送 | `trivEquiv_restrict`（在庫、`LocalMetric.lean`） |
| 遷移単元で座標が動く | `trivEquiv_transUnit`（本ファイル、`trivValue_transUnit` の切断版） |

## ★残っている段（明示）

★★本ファイルは**判定基準だけ**である。段 E3a-3 を閉じるには、
`§9-822`（分母を払う）が出す `a_j` がこの判定を通ることを示す段が要る
——そこで `§9-826`（重なりの一致は冪で強制できる）が効く。
★★★貼り合わせ自体は `§9-827`（`exists_unique_glue_top`）が持っている。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Opposite MonoidalCategory ABC3.Found.Arakelov

variable {X : Scheme.{0}}

/-! ## ★遷移単元は `trivEquiv` を動かす -/

/-- ★★**遷移単元は `trivEquiv` を動かす**（`trivValue_transUnit` の「`V` 上の切断」版）。

★在庫の `trivValue_transUnit` は**大域切断**についての形だったが、
貼り合わせでは `V` 上の切断（制限したもの）を比べるので、こちらの形が要る。 -/
theorem trivEquiv_transUnit (F : X.PresheafOfModules) (V : X.Opens)
    (e e' : (restrictPresheafFunctor X V).obj F ≅ 𝟙_ (PresheafModulesOn X V))
    (x : (F.obj (op V) : Type)) :
    trivEquiv F V e' x = trivEquiv F V e x * transUnit F V e e' := by
  have h := linearEquiv_self_apply
    ((trivEquiv F V e).symm.trans (trivEquiv F V e')) (trivEquiv F V e x)
  show _ = _ * ((trivEquiv F V e).symm.trans (trivEquiv F V e')) 1
  simpa using h

/-! ## ★★★★テンソル冪の自明化は制限・遷移と両立する -/

set_option backward.isDefEq.respectTransparency false in
/-- ★★★**`tensorPowTriv` は制限と可換である**。

★各段が `rfl`（`trivialOfLe` は `tensorTriv` と可換）なので、再帰で出る。 -/
theorem trivialOfLe_tensorPowTriv {M : X.PresheafOfModules} {V W : X.Opens} (hWV : W ≤ V)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V)) : ∀ n : ℕ,
    trivialOfLe (presheafTensorPow M n) hWV (tensorPowTriv e n)
      = tensorPowTriv (trivialOfLe M hWV e) n
  | 0 => rfl
  | n + 1 => by
      show trivialOfLe (M ⊗ presheafTensorPow M n) hWV (tensorTriv e (tensorPowTriv e n)) = _
      show tensorTriv (trivialOfLe M hWV e)
        (trivialOfLe (presheafTensorPow M n) hWV (tensorPowTriv e n)) = _
      rw [trivialOfLe_tensorPowTriv hWV e n]
      rfl

/-- ★★★★★★★**`M^{⊗n}` の遷移単元は `M` の遷移単元の `n` 乗である**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★★これが「`n` 冪へ移る」ことの**幾何的な意味**である
——遷移函数が `n` 乗されるから、`n` を上げれば分母の極が飲み込める。

★★★機構は `transUnit_tensorTriv`（在庫——遷移単元はテンソル積で掛け算になる）の再帰と、
基点の `transUnit_self`（同じ自明化どうしの遷移単元は `1`）だけである。 -/
theorem transUnit_tensorPowTriv {M : X.PresheafOfModules} {V : X.Opens}
    (e e' : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V)) : ∀ n : ℕ,
    transUnit (presheafTensorPow M n) V (tensorPowTriv e n) (tensorPowTriv e' n)
      = (transUnit M V e e') ^ n
  | 0 => by
      show transUnit (𝟙_ X.PresheafOfModules) V restrictPresheafUnit.symm
        restrictPresheafUnit.symm = _
      rw [pow_zero, transUnit_self]
  | n + 1 => by
      show transUnit (M ⊗ presheafTensorPow M n) V
        (tensorTriv e (tensorPowTriv e n)) (tensorTriv e' (tensorPowTriv e' n)) = _
      rw [transUnit_tensorTriv, transUnit_tensorPowTriv e e' n, pow_succ]
      ring

/-! ## ★★★★★関数から切断へ -/

/-- ★★★★**チャートの関数 `a ∈ Γ(X, V)` が定める `M^{⊗N}` の切断**。

★段 E3a-3 で「分母を払って得た関数」を切断に戻すのがこれである。 -/
noncomputable def secOfFun (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V)) (N : ℕ)
    (a : (Γ(X, V) : Type)) : ((presheafTensorPow M N).obj (op V) : Type) :=
  (trivEquiv (presheafTensorPow M N) V (tensorPowTriv e N)).symm a

/-- ★**座標を読むと元の関数に戻る**。 -/
theorem trivEquiv_secOfFun (M : X.PresheafOfModules) (V : X.Opens)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V)) (N : ℕ)
    (a : (Γ(X, V) : Type)) :
    trivEquiv (presheafTensorPow M N) V (tensorPowTriv e N) (secOfFun M V e N a) = a :=
  LinearEquiv.apply_symm_apply _ a

/-- ★★★**制限しても「関数から作った切断」のままである**——関数の側も制限されるだけ。 -/
theorem trivEquiv_res_secOfFun {M : X.PresheafOfModules} {V W : X.Opens} (hWV : W ≤ V)
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V)) (N : ℕ)
    (a : (Γ(X, V) : Type)) :
    trivEquiv (presheafTensorPow M N) W (tensorPowTriv (trivialOfLe M hWV e) N)
        ((presheafTensorPow M N).map (homOfLE hWV).op (secOfFun M V e N a))
      = X.presheaf.map (homOfLE hWV).op a := by
  rw [secOfFun, ← trivialOfLe_tensorPowTriv hWV e N,
    trivEquiv_restrict (presheafTensorPow M N) hWV (tensorPowTriv e N),
    LinearEquiv.apply_symm_apply]

/-! ## ★★★★★★★★判定基準 -/

/-- ★★★★★★★★**重なりの一致判定** —— 段 E3a-3 の判定基準。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

    `secOfFun e N a` と `secOfFun e' N a'` が `V ⊓ V'` で一致
      ⟺  `a|_{V⊓V'} · u^N = a'|_{V⊓V'}`

ここで `u = transUnit M (V ⊓ V') (e|_{V⊓V'}) (e'|_{V⊓V'})` は遷移単元である。

★★**遷移単元がちょうど `N` 乗で効く**のが要点である
——`M^{⊗N}` の遷移函数は `M` の遷移函数の `N` 乗だからである（`transUnit_tensorPowTriv`）。

★★★これで「切断が貼り合うか」という幾何の問いが、
**`Γ(X, V ⊓ V')` の中の等式**という代数の問いに落ちる。 -/
theorem res_secOfFun_eq_iff {M : X.PresheafOfModules} {V V' : X.Opens}
    (e : (restrictPresheafFunctor X V).obj M ≅ 𝟙_ (PresheafModulesOn X V))
    (e' : (restrictPresheafFunctor X V').obj M ≅ 𝟙_ (PresheafModulesOn X V'))
    (N : ℕ) (a : (Γ(X, V) : Type)) (a' : (Γ(X, V') : Type)) :
    (presheafTensorPow M N).map (homOfLE (inf_le_left : V ⊓ V' ≤ V)).op
        (secOfFun M V e N a)
      = (presheafTensorPow M N).map (homOfLE (inf_le_right : V ⊓ V' ≤ V')).op
        (secOfFun M V' e' N a')
    ↔ X.presheaf.map (homOfLE (inf_le_left : V ⊓ V' ≤ V)).op a
        * (transUnit M (V ⊓ V')
            (trivialOfLe M (inf_le_left : V ⊓ V' ≤ V) e)
            (trivialOfLe M (inf_le_right : V ⊓ V' ≤ V') e')) ^ N
      = X.presheaf.map (homOfLE (inf_le_right : V ⊓ V' ≤ V')).op a' := by
  have key : ∀ x : ((presheafTensorPow M N).obj (op (V ⊓ V')) : Type),
      trivEquiv (presheafTensorPow M N) (V ⊓ V')
        (tensorPowTriv (trivialOfLe M (inf_le_right : V ⊓ V' ≤ V') e') N) x
      = trivEquiv (presheafTensorPow M N) (V ⊓ V')
          (tensorPowTriv (trivialOfLe M (inf_le_left : V ⊓ V' ≤ V) e) N) x
        * (transUnit M (V ⊓ V')
            (trivialOfLe M (inf_le_left : V ⊓ V' ≤ V) e)
            (trivialOfLe M (inf_le_right : V ⊓ V' ≤ V') e')) ^ N := by
    intro x
    rw [trivEquiv_transUnit _ _
      (tensorPowTriv (trivialOfLe M (inf_le_left : V ⊓ V' ≤ V) e) N),
      transUnit_tensorPowTriv]
  constructor
  · intro h
    have h2 := congrArg (fun x => trivEquiv (presheafTensorPow M N) (V ⊓ V')
      (tensorPowTriv (trivialOfLe M (inf_le_right : V ⊓ V' ≤ V') e') N) x) h
    simp only at h2
    rw [trivEquiv_res_secOfFun (inf_le_right : V ⊓ V' ≤ V') e' N a', key,
      trivEquiv_res_secOfFun (inf_le_left : V ⊓ V' ≤ V) e N a] at h2
    exact h2
  · intro h
    refine (trivEquiv (presheafTensorPow M N) (V ⊓ V')
      (tensorPowTriv (trivialOfLe M (inf_le_right : V ⊓ V' ≤ V') e') N)).injective ?_
    rw [trivEquiv_res_secOfFun (inf_le_right : V ⊓ V' ≤ V') e' N a', key,
      trivEquiv_res_secOfFun (inf_le_left : V ⊓ V' ≤ V) e N a]
    exact h

/-! ## ★出典の紐付け(`.src`) -/

def secOfFun.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(チャートの関数が定める M^{⊗N} の切断)",
    sectionId := "genell-prop-1-4" }

def transUnit_tensorPowTriv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(M^{⊗n} の遷移単元は M の遷移単元の n 乗)",
    sectionId := "genell-prop-1-4" }

def trivialOfLe_tensorPowTriv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(tensorPowTriv は制限と可換)",
    sectionId := "genell-prop-1-4" }

def res_secOfFun_eq_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(重なりの一致判定——遷移単元の N 乗)",
    sectionId := "genell-prop-1-4" }

def res_secOfFun_eq_iff.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "transUnit_tensorTriv(遷移単元はテンソル積で掛け算になる)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.transUnit_tensorTriv") 3,
    .citation "[ABC3]" "trivEquiv_restrict(座標は制限と可換、LocalMetric.lean)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.trivEquiv_restrict") 3,
    .citation "[Stacks]" "Lemma 01PW(の貼り合わせの段)"
      (.absent "mathlib に ample は無い(2026-08-28 実測)") 7,
    .implicitStep
      ("★★**遷移単元がちょうど N 乗で効く**のが「n 冪へ移る」ことの幾何的な意味である" ++
       "——M^{⊗N} の遷移函数は M の遷移函数の N 乗だから、N を上げれば分母の極が飲み込める。" ++
       "★原文はこれを「[some positive tensor power of]」と 1 語で畳んでいる") 8,
    .implicitStep
      ("★★★本ファイルは**判定基準だけ**である。段 E3a-3 を閉じるには、" ++
       "§9-822(分母を払う)が出す a_j がこの判定を通ることを示す段が要る" ++
       "——そこで §9-826(重なりの一致は冪で強制できる)が効く。" ++
       "貼り合わせ自体は §9-827(exists_unique_glue_top)が持っている") 8 ]

end ABC3.Found.GenEll
